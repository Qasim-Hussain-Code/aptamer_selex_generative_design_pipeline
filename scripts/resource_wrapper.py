"""Sample process-tree RSS and all project files; terminate before resource exhaustion."""
import argparse
import os
import shlex
import subprocess
import time
import shutil
from common import ROOT, append, disk_bytes, failure, guard, now


def main():
    import psutil
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", required=True)
    p.add_argument("command", nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command[:1] == ["--"] else a.command
    guard()
    start, wall = now(), time.monotonic()
    before = disk_bytes()
    peak_disk, peak_rss, minimum_free = before, 0, shutil.disk_usage(ROOT).free
    proc = subprocess.Popen(command, cwd=ROOT)
    root_proc = psutil.Process(proc.pid)
    note = "RSS sums sampled process tree at 0.2 s; disk sampled at 1 s; short peaks may be missed"
    last_disk = 0
    stopped = False
    while proc.poll() is None:
        try:
            family = [root_proc] + root_proc.children(recursive=True)
            rss = sum(q.memory_info().rss for q in family if q.is_running())
            peak_rss = max(peak_rss, rss)
            if time.monotonic() - last_disk >= 1:
                peak_disk = max(peak_disk, disk_bytes())
                minimum_free = min(minimum_free, shutil.disk_usage(ROOT).free)
                last_disk = time.monotonic()
            if rss > int(os.environ.get("RAM_LIMIT_BYTES", 14_000_000_000)) or peak_disk > int(os.environ.get("DISK_BUDGET_BYTES", 13_000_000_000)) or minimum_free < int(os.environ.get("RESERVE_BYTES", 1_000_000_000)):
                stopped = True
                for q in reversed(family):
                    try:
                        q.terminate()
                    except psutil.Error:
                        pass
                note += "; terminated for resource ceiling"
                break
        except psutil.Error:
            pass
        time.sleep(0.2)
    try:
        status = proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        proc.kill()
        status = proc.wait()
    after = disk_bytes()
    append(ROOT / "logs/resource_usage.tsv", dict(stage=a.stage, command=shlex.join(command),
           start_timestamp=start, end_timestamp=now(), elapsed_seconds=f"{time.monotonic()-wall:.3f}",
           exit_status=status, peak_rss_bytes=peak_rss, disk_bytes_before=before,
           peak_observed_disk_bytes=max(peak_disk, after), disk_bytes_after=after,
           threads=os.environ.get("THREADS", "1"), execution_host_class="local_cpu",
           notes=note, minimum_free_bytes=minimum_free))
    if status or stopped:
        failure(a.stage, shlex.join(command), RuntimeError(f"exit={status}; {note}"))
    raise SystemExit(status if status else (75 if stopped else 0))


if __name__ == "__main__":
    main()
