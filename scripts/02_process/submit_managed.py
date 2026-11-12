#!/usr/bin/env python
"""submit_managed.py -- one SLURM job per group, driven by InSARHub's manager.

Replaces the job-array lanes. An array bundles many groups under one job id and
one resource shape, which makes the campaign hard to steer: a single slow group
holds a task slot, `scontrol` can only move the whole array, and the `all`
partition's QOS counts array tasks as submissions so the array itself has to be
strided to stay legal.

InSARHub already solves this, in utils/slurm_manager.py: a lightweight manager
job keeps at most ``max_concurrent`` CHILD jobs alive, each its own sbatch
script with its own job id, and refills a slot the moment one finishes. It
handles the parts that are easy to get wrong -- it treats
QOSMaxSubmitJobPerUserLimit as backpressure and retries rather than dying, and
it waits out NFS attribute lag before believing a child produced no marker,
falling back to sacct. None of that is reimplemented here.

LANES. The managers all run on smi_all (30-day walltime, no QOS), and their
children are split by what the partitions impose:

    uncapped   children on smi_all,smi+coe   -- QoS=N/A on both, so no ceiling
    capped     children on all,coe_all,coe_cpu -- all three carry
               QoS=jobmanagement, MaxSubmitJobsPerUser=20

The split is not cosmetic. sbatch applies the STRICTEST QOS across a job's
partition list, so a single manager whose children named every partition would
put the uncapped work under the 20-job ceiling too. Two managers, each with its
own partition list, is what lets the capped family sit at exactly its 20 while
the uncapped family runs as wide as the cluster allows.

The manager itself costs one idle core and is not subject to the ceiling, since
smi_all imposes no QOS -- which is why the whole campaign can be steered from a
partition that has no submit limit.

Group lists per lane must be disjoint: two jobs in the same group directory
would race on cslc/ and pairs/.

Usage:
    python submit_managed.py --root <wx root> --list plan/v4_groups.txt \
        [--exclude GID ...] [--uncapped 96] [--dry-run]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/home/C838053462/dev/InSARHub/src")
from insarhub.utils.slurm_manager import (            # noqa: E402
    build_sliding_window_manager, qos_max_submit_jobs,
)

CHILD = "/home/C838053462/dev/extreme_wx/slurm/run_one_group.sh"
ENV = "/home/C838053462/miniforge3/envs/ifgqual/lib"

# A group is one burst: two acquisitions to download and geocode and one
# interferogram, measured at ~20 min once the date list is pruned. 8 h is the
# margin for a slow ASF transfer, and is legal on every partition here (`all`
# and smi+coe cap at 2 days, smi_all at 30).
# 45 min against a measured 7 min 9 s, and the margin is deliberately small.
#
# This account's jobs sit at fair-share priority 1, the cluster floor, so they
# never start by priority -- they start by BACKFILL, and backfill only runs a
# job if it will finish before the higher-priority reservation it is filling in
# front of. An 8-hour request is too long to fit almost any gap, which is why
# 105 of our jobs sat queued while engineering-cpu1-4 held 60 idle cores. The
# walltime is a scheduling parameter here, not a safety margin: asking for more
# than the work needs buys nothing and costs every backfill opportunity.
CHILD_TIME = "00:45:00"
# 2 cores, not 4. A pruned group is 1 burst x 2 dates, so the cslc stage has
# exactly two parallel jobs and the ifg stage one -- `--worker 4` can never use
# more than two of them. The 4 was sized for the main campaign's 4-to-8-date
# groups and carried over by habit. Halving it doubles how many jobs fit a node
# (smith-gpu1 went from 15 to ~31) and makes each one small enough to backfill
# into the 12-16 core gaps on the engineering nodes, which the scheduler was
# otherwise holding for higher-priority work while our 96 queued jobs waited.
CHILD_CPUS = 2
# 8G against a measured 288 MB peak. Memory is not the binding constraint here,
# but it is a scheduling one: 48G is what held the main campaign to 6 running
# tasks on a 128-core partition, and anything that does not fit a node's spare
# memory waits just as surely as one that does not fit its spare cores.
CHILD_MEM = "8G"

LANES = {
    # name:      (child partitions, default max_concurrent)
    "uncapped": ("smi_all,smi+coe", 96),
    "capped":   ("all,coe_all,coe_cpu", None),   # None -> read the QOS ceiling
}
MANAGER_PARTITION = "smi_all"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--list", required=True)
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="group ids to leave out (e.g. one already running)")
    ap.add_argument("--uncapped", type=int, default=LANES["uncapped"][1],
                    help="concurrent child jobs on smi_all,smi+coe")
    ap.add_argument("--tag", default="v4")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    lp = Path(a.list)
    gids = [g for g in (lp if lp.is_absolute() else a.root / lp).read_text().split() if g]
    excl = {e.upper() for e in a.exclude}
    gids = [g for g in gids if g.upper() not in excl]
    if excl:
        print(f"excluded {len(excl)} group(s): {', '.join(sorted(excl))}")

    cap, qos = qos_max_submit_jobs("all")
    # Leave one slot of the ceiling free. The manager retries a refused
    # submission rather than failing, so running right at 20 only costs it
    # polling -- but the ceiling is USER-WIDE, and sitting on all of it means
    # any other job this account submits is refused for as long as the
    # campaign runs.
    capped_n = max(1, (cap or 20) - 1)
    print(f"QOS on `all`: {qos}, MaxSubmitJobsPerUser={cap} -> capped lane at {capped_n}")

    slots = {"uncapped": a.uncapped, "capped": capped_n}
    total = sum(slots.values())
    # Deal round-robin weighted by each lane's slot count, so both lanes finish
    # at about the same time and neither is left holding the tail.
    import itertools
    wheel = list(itertools.chain.from_iterable([k] * slots[k] for k in LANES))
    lanes: dict[str, list[str]] = {k: [] for k in LANES}
    for i, g in enumerate(gids):
        lanes[wheel[i % len(wheel)]].append(g)

    seen: set[str] = set()
    for k, v in lanes.items():
        assert not (set(v) & seen), f"lane {k} overlaps another"
        seen |= set(v)
    assert len(seen) == len(gids)

    hpc = a.root / "_campaign" / f"managed_{a.tag}"
    submitted = []
    for name, (partitions, _) in LANES.items():
        gl = lanes[name]
        if not gl:
            continue
        d = hpc / name
        (d / "logs").mkdir(parents=True, exist_ok=True)
        (d / "sbatch").mkdir(parents=True, exist_ok=True)
        status = d / "status"
        status.mkdir(exist_ok=True)

        cmds = [f"bash {CHILD} {g}" for g in gl]
        write_status_fn = (
            f'STATUS_DIR="{status}"\n'
            'write_status() {\n'
            '    rm -f "$STATUS_DIR/.succeeded" "$STATUS_DIR/.failed"\n'
            '    if [[ "$1" == "SUCCEEDED" ]]; then touch "$STATUS_DIR/.succeeded";\n'
            '    else touch "$STATUS_DIR/.failed"; fi\n'
            '    echo "$1" > "$STATUS_DIR/status.txt"\n'
            '}'
        )
        mgr = build_sliding_window_manager(
            job_name_base=f"wx{a.tag}_{name}"[:32],
            commands=cmds,
            log_dir=d / "logs",
            sbatch_dir=d / "sbatch",
            max_concurrent=slots[name],
            slurm_kwargs=dict(partition=partitions, time=CHILD_TIME,
                              cpus_per_task=CHILD_CPUS, mem=CHILD_MEM,
                              nodes=1, ntasks=1),
            env_lines=[
                f'export LD_LIBRARY_PATH={ENV}${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}',
                f'export WX_ROOT={a.root.resolve()}',
                'export WX_SLURM_DIR=/home/C838053462/dev/extreme_wx/slurm',
            ],
            write_status_fn=write_status_fn,
            # The manager outlives every child it supervises, so it takes the
            # 30-day partition rather than the 2-day ones the work runs on.
            manager_partition=MANAGER_PARTITION,
            file_prefix=name, label=f"{a.tag}/{name}",
        )
        print(f"{name:9s} {len(gl):>5,} groups  {slots[name]:>3} concurrent  "
              f"children -> {partitions}")
        if a.dry_run:
            print(f"          manager script: {mgr}")
            continue
        r = subprocess.run(["sbatch", str(mgr)], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  sbatch FAILED: {r.stderr.strip()}", file=sys.stderr)
            continue
        jid = r.stdout.strip().split()[-1]
        submitted.append((name, jid, len(gl)))
        print(f"          manager job {jid}")

    if submitted and not a.dry_run:
        print()
        print(f"{len(gids):,} groups, {total} concurrent child jobs at peak")
        for name, jid, n in submitted:
            print(f"  manager {jid}  {name}  {n:,} groups")


if __name__ == "__main__":
    main()
