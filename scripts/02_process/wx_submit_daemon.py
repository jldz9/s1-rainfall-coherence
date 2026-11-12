#!/usr/bin/env python
"""wx_submit_daemon.py -- sliding-window submitter: one SLURM job per group.

Same design as InSARHub's utils/slurm_manager.py -- keep N child jobs alive,
refill a slot the moment one leaves -- but driven from a process on the login
node instead of from a manager job inside the cluster.

The manager-job version is the better design in general and it is what was
tried first: jobs 607401 and 607402 are two of them, and both are still
queued, because smi_all's 128 cores are held by another user's 64-core job and
a manager cannot submit anything until it is itself scheduled. A manager that
needs a node to hand out nodes stalls the whole campaign behind one core. This
one needs nothing but squeue and sbatch.

Two differences from the library version, both deliberate:

  * It does not stop on a failed child. The library gates its refill on
    FAIL_COUNT == 0, which is right for a staged pipeline where a failed cslc
    makes the following ifg meaningless, and wrong for 1,464 independent
    groups of which a handful carry known data defects -- one of those would
    strand the rest.
  * Its pool targets are re-read from disk every poll, so concurrency can be
    retuned while it runs, the way `scontrol update ArrayTaskThrottle` retunes
    a live array. Raising a target takes effect on the next poll; lowering one
    never cancels anything, it just stops refilling until the pool drains to
    the new size.

POOLS. One per QOS family, because sbatch applies the strictest QOS across a
job's partition list -- a single pool naming every partition would put the
uncapped work under the capped family's ceiling:

    u   smi_all,smi+coe        QoS=N/A on both, so no submit ceiling
    c   all,coe_all,coe_cpu    all three carry QoS=jobmanagement,
                               MaxSubmitJobsPerUser=20, counted USER-WIDE

The pool is encoded in each job's name (``u_<group>`` / ``c_<group>``), so the
daemon's whole state is recoverable from squeue and it can be restarted
without losing track of what is in flight.

No group list is held in memory as the source of truth: a group is work if it
has no .done, and run_one_group.sh takes a per-group flock, so a second
submitter is wasteful rather than destructive.

Usage:
    python wx_submit_daemon.py --root <wx root> --list plan/v4_groups.txt \
        [--uncapped 96] [--capped 19] [--poll 20] [--once]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

CHILD = "/home/C838053462/dev/extreme_wx/slurm/run_one_group.sh"
# smi+coe was moved out of the uncapped pool on 2026-10-08. sbatch applies the
# strictest QoS across a job's partition list, and smi+coe acquired
# QoS=jobmanagement (MaxSubmitJobsPerUser=20) part way through the campaign, so
# every job named for both partitions inherited that ceiling. The daemon then
# hit QOSMaxSubmitJobPerUserLimit two submissions into each pass and broke out
# of its loop, leaving a third of smi_all idle while the pool drained.
#
# smi_all carries QoS=N/A and 128 cores, so on its own it has no submit ceiling.
POOLS = {
    "u": "smi_all",
    "c": "all,coe_all,coe_cpu,smi+coe",
}
# sbatch refuses a submission over the QOS ceiling with QOSMaxSubmitJobPerUserLimit.
# That is backpressure, not an error: it clears as earlier children finish, so it
# is logged once and retried on the next poll rather than treated as a failure.
BACKPRESSURE = ("QOSMaxSubmitJobPerUserLimit", "AssocMaxSubmitJobLimit",
                "job submit limit")


def squeue_mine() -> dict[str, str]:
    """{job_id: job_name} for this user's active jobs, or {} if squeue fails.

    An empty dict on failure is deliberately indistinguishable from an empty
    queue here: the daemon would then submit up to its targets again. That is
    safe only because run_one_group.sh locks its group -- a duplicate exits
    immediately instead of corrupting the group it landed on.
    """
    try:
        r = subprocess.run(["squeue", "--me", "--noheader", "--format=%i %j"],
                           capture_output=True, text=True, timeout=30)
        out = {}
        for line in r.stdout.splitlines():
            p = line.split(None, 1)
            if len(p) == 2:
                out[p[0]] = p[1].strip()
        return out
    except Exception as e:                                       # noqa: BLE001
        print(f"  squeue failed ({type(e).__name__}: {e}); holding this poll", flush=True)
        return {}


def submit(gid: str, pool: str, root: Path, log_dir: Path,
           walltime: str = "00:59:00") -> tuple[str | None, str]:
    # Resources come from the command line, not from run_one_group.sh's own
    # #SBATCH header: that file is being executed by jobs already on nodes, and
    # bash reads a script incrementally, so editing it would make them resume at
    # a stale byte offset (README trap 9). sbatch flags override the header and
    # touch nothing that is running.
    #
    # 45 min, not 8 h, because priority here is PriorityWeightAge alone
    # (FairShare, JobSize, QOS and Partition are all weighted 0) and everything
    # we start, we start by backfill. bf_resolution=3600 means backfill plans in
    # whole hours, so 45 min and 59 min cost the same and 8 h costs eight times
    # as much. 59 minutes rather than 45 for exactly that reason: both round to
    # the same whole hour and cost the same, and at 45 the measured mean group
    # time of 34.7 min left so little margin that 81 of the first 150 finished
    # jobs were killed at the wall.
    #
    # 1 core and 2 GB: measured on an identical group at seven resource
    # levels, memory made no difference at all (15:22 at 2 GB against 15:25 at 4,
    # 9:43 at 6 GB against 9:34 at 8) and cores scaled sub-linearly (7:50 on 4,
    # 9:43 on 2, 15:22 on 1), so one core per job finishes the same work in
    # 2.04x the throughput per core that four do, and 2 GB never came close to
    # being exhausted.
    cmd = ["sbatch",
           f"--partition={POOLS[pool]}",
           f"--job-name={pool}_{gid}",
           f"--time={walltime}", "--cpus-per-task=1", "--mem=2G",
           f"--output={log_dir}/{gid}_%j.out",
           f"--error={log_dir}/{gid}_%j.err",
           f"--export=ALL,WX_ROOT={root},"
           f"WX_SLURM_DIR=/home/C838053462/dev/extreme_wx/slurm,"
           f"LD_LIBRARY_PATH=/home/C838053462/miniforge3/envs/ifgqual/lib",
           CHILD, gid]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except Exception as e:                                       # noqa: BLE001
        return None, f"{type(e).__name__}: {e}"
    if r.returncode == 0:
        return r.stdout.strip().split()[-1], ""
    return None, (r.stderr or r.stdout).strip()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--list", required=True)
    ap.add_argument("--uncapped", type=int, default=96,
                    help="in-flight jobs on smi_all,smi+coe (384 CPUs / 4 per job)")
    ap.add_argument("--capped", type=int, default=19,
                    help="in-flight jobs under QoS=jobmanagement (ceiling is 20)")
    ap.add_argument("--poll", type=int, default=20)
    ap.add_argument("--walltime", default="00:59:00",
                    help="worker wall time. 59 min was chosen when backfill was "
                         "contested; with smi_all effectively ours, a longer "
                         "limit costs little and removes the timeout cliff that "
                         "killed 81 of the first 150 jobs at 45 min")
    ap.add_argument("--once", action="store_true", help="one pass, then exit")
    a = ap.parse_args()

    root = a.root.resolve()
    lp = Path(a.list)
    gids = [g for g in (lp if lp.is_absolute() else root / lp).read_text().split() if g]

    log_dir = root / "_campaign" / "managed_v4" / "jobs"
    log_dir.mkdir(parents=True, exist_ok=True)
    # The CLI flags seed the file on every run. Honouring a stale file instead
    # silently ignored `--uncapped 836` and submitted 96, because an earlier
    # invocation had written 96 into it. The file's job is live RETUNING of a
    # running daemon, not remembering what the last one was told.
    targets_file = root / "_campaign" / "managed_v4" / "targets.json"
    targets_file.write_text(json.dumps({"u": a.uncapped, "c": a.capped}, indent=2))

    print(f"[{time.strftime('%F %T')}] daemon start: {len(gids):,} groups in list", flush=True)
    print(f"  targets from {targets_file} (edit it to retune while running)", flush=True)

    announced_backpressure: set[str] = set()
    pass_n = 0
    while True:
        pass_n += 1
        try:
            targets = json.loads(targets_file.read_text())
        except Exception:                                        # noqa: BLE001
            targets = {"u": a.uncapped, "c": a.capped}

        active = squeue_mine()
        inflight = {p: set() for p in POOLS}
        for name in active.values():
            if len(name) > 2 and name[1] == "_" and name[0] in POOLS:
                inflight[name[0]].add(name[2:])
        # Only this run's own groups count: another submitter's jobs are real
        # but they are not progress against THIS list, and counting them made
        # `done` go negative.
        gset = set(gids)
        busy = {g for s in inflight.values() for g in s} & gset

        todo = [g for g in gids
                if g not in busy and not (root / "proc" / "burst" / g / ".done").exists()]
        done = len(gids) - len(todo) - len(busy)

        subs = 0
        for pool in POOLS:
            want = max(0, int(targets.get(pool, 0)) - len(inflight[pool]))
            while want > 0 and todo:
                gid = todo.pop(0)
                jid, err = submit(gid, pool, root, log_dir, a.walltime)
                if jid:
                    subs += 1
                    want -= 1
                    inflight[pool].add(gid)
                    continue
                if any(b in err for b in BACKPRESSURE):
                    if pool not in announced_backpressure:
                        print(f"  pool {pool}: at the queue submit limit; "
                              f"will retry each poll", flush=True)
                        announced_backpressure.add(pool)
                    todo.insert(0, gid)
                    break
                print(f"  sbatch FAILED {gid} ({pool}): {err}", flush=True)
                # Not retried this pass; it stays in the list and comes back on
                # the next poll, so a transient scheduler error self-heals and a
                # permanent one is visible in the log every poll instead of once.
                break

        print(f"[{time.strftime('%F %T')}] pass {pass_n}: "
              f"done {done:,}/{len(gids):,} | running+queued u={len(inflight['u'])} "
              f"c={len(inflight['c'])} | submitted {subs} | remaining {len(todo):,}",
              flush=True)

        if not todo and not busy:
            print(f"[{time.strftime('%F %T')}] all groups accounted for; daemon exiting",
                  flush=True)
            break
        if a.once:
            break
        time.sleep(a.poll)


if __name__ == "__main__":
    sys.exit(main())
