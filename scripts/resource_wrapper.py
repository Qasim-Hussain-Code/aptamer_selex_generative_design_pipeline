"""Sample process-tree RSS and all project files; terminate before resource exhaustion."""
import argparse
import os
import shlex
import subprocess
import time
import shutil
from common import ROOT, append, disk_bytes, failure, guard, now


def main():
    try:
        import psutil
    except ImportError:
        psutil=None
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", required=True)
    p.add_argument("command", nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command[:1] == ["--"] else a.command
    guard()
    start, wall = now(), time.monotonic()
    before = disk_bytes()
    # The environment is immutable outside installation. Its measured bytes still count,
    # but repeated stat calls on thousands of package files would dominate short stages.
    frozen_env = ROOT / ".venv"
    frozen_bytes = disk_bytes(frozen_env) if frozen_env.exists() and a.stage != "install" else 0
    def observed_disk():
        if not frozen_bytes:
            return disk_bytes()
        total=frozen_bytes
        for parent, directories, files in os.walk(ROOT):
            if os.path.abspath(parent)==str(ROOT):
                directories[:]=[d for d in directories if d!=".venv"]
            for name in files:
                try:
                    total+=(__import__("pathlib").Path(parent)/name).stat().st_size
                except FileNotFoundError:
                    pass
        return total
    peak_disk, peak_rss, minimum_free = before, 0, shutil.disk_usage(ROOT).free
    proc = subprocess.Popen(command, cwd=ROOT)
    root_proc = psutil.Process(proc.pid) if psutil else None
    note = "RSS sums sampled process tree at 0.2 s; disk sampled at 1 s; short peaks may be missed"
    last_disk = 0
    stopped = False
    rss_errors = set()
    while proc.poll() is None:
        if psutil is None:
            if os.name == "nt":
                import ctypes
                from ctypes import wintypes
                class MemoryCounters(ctypes.Structure):
                    _fields_=[("cb",wintypes.DWORD),("PageFaultCount",wintypes.DWORD)]+[(n,ctypes.c_size_t) for n in ["PeakWorkingSetSize","WorkingSetSize","QuotaPeakPagedPoolUsage","QuotaPagedPoolUsage","QuotaPeakNonPagedPoolUsage","QuotaNonPagedPoolUsage","PagefileUsage","PeakPagefileUsage"]]
                counter=MemoryCounters();counter.cb=ctypes.sizeof(counter)
                ok=ctypes.windll.psapi.GetProcessMemoryInfo(wintypes.HANDLE(int(proc._handle)),ctypes.byref(counter),counter.cb)
                if ok:peak_rss=max(peak_rss,counter.PeakWorkingSetSize)
            elif os.path.exists(f"/proc/{proc.pid}/status"):
                with open(f"/proc/{proc.pid}/status") as f:
                    for line in f:
                        if line.startswith("VmHWM:"):peak_rss=max(peak_rss,int(line.split()[1])*1024)
            peak_disk=max(peak_disk,observed_disk())
            minimum_free=min(minimum_free,shutil.disk_usage(ROOT).free)
            if peak_rss>int(os.environ.get("RAM_LIMIT_BYTES",14_000_000_000)) or peak_disk>int(os.environ.get("DISK_BUDGET_BYTES",13_000_000_000)) or minimum_free<int(os.environ.get("RESERVE_BYTES",1_000_000_000)):
                proc.terminate();stopped=True;break
            time.sleep(0.2)
            continue
        try:
            family = [root_proc]
            try:
                family += root_proc.children(recursive=True)
            except psutil.Error as e:
                rss_errors.add(f"children:{type(e).__name__}")
            rss = 0
            for q in family:
                try:
                    rss += q.memory_info().rss
                except psutil.Error as e:
                    rss_errors.add(f"memory:{type(e).__name__}")
            peak_rss = max(peak_rss, rss)
            if time.monotonic() - last_disk >= 1:
                peak_disk = max(peak_disk, observed_disk())
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
    if psutil is None:
        note += "; bootstrap native process peak RSS; child-tree RSS unavailable until psutil installed"
    if rss_errors:
        note += "; RSS sampling errors=" + ",".join(sorted(rss_errors))
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
