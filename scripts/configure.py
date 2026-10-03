import argparse
import json
import os
import platform
import shlex
import shutil
from pathlib import Path
from common import ROOT, atomic_text, disk_bytes, table, now


def effective_budget(requested, free, reserve=1_000_000_000):
    if free <= reserve:
        raise ValueError("Free space does not exceed the external reserve")
    return min(int(requested), 13_000_000_000, free - reserve)


def main():
    p = argparse.ArgumentParser(description="Measure resources and write shell configuration")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--ram", type=float, default=16)
    p.add_argument("--disk", type=float, default=13)
    p.add_argument("--gpu-mode", choices=["none", "hosted", "local"], default="none")
    p.add_argument("--yes", action="store_true")
    a = p.parse_args()
    if a.threads < 1 or a.ram <= 2 or a.disk <= 0:
        p.error("threads must be positive, RAM above 2 GB, disk positive")
    free = shutil.disk_usage(ROOT).free
    budget = effective_budget(a.disk * 1e9, free)
    if disk_bytes() >= budget:
        raise RuntimeError("Existing project exceeds effective budget")
    values = dict(THREADS=a.threads, RAM_LIMIT_BYTES=int(min(a.ram - 2, 14) * 1e9),
                  DISK_BUDGET_BYTES=budget, RESERVE_BYTES=1_000_000_000,
                  GPU_MODE=a.gpu_mode, KEEP_CHECKPOINTS="false")
    atomic_text(ROOT / "project.conf", "\n".join(f"export {k}={shlex.quote(str(v))}" for k, v in values.items()) + "\n")
    table(ROOT / "results/resource_configuration.tsv", [dict(timestamp=now(), free_bytes=free,
          requested_disk_bytes=int(a.disk * 1e9), effective_budget_bytes=budget,
          reserve_bytes=1_000_000_000, project_bytes=disk_bytes(), ram_limit_bytes=values["RAM_LIMIT_BYTES"], threads=a.threads)])
    print(f"Configured budget {budget} bytes, free {free} bytes, reserve 1000000000 bytes")


if __name__ == "__main__":
    main()
