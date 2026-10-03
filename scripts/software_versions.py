import importlib.metadata as im
import json
import platform
import subprocess
import sys
import os
import psutil
from common import ROOT, table, now, atomic_text


def main():
    records = []
    for line in (ROOT / "config/requirements.txt").read_text().splitlines():
        name, pinned = line.split("==")
        dist = im.distribution(name)
        licence = dist.metadata.get("License-Expression") or dist.metadata.get("License") or "unclear"
        if not licence or licence == "UNKNOWN":
            licence = "; ".join(s for s in dist.metadata.get_all("Classifier", []) if "License" in s) or "unclear"
        records.append(dict(tool=name, role="local core", version=dist.version,
             source=f"https://pypi.org/project/{name}/{dist.version}/", source_commit="", licence=licence,
             licence_url_or_source="installed distribution METADATA and licence files", licence_checked_date="2026-10-03",
             execution="local", redistributable="yes" if licence != "unclear" else "unclear", notes="Installed separately; no dependency source vendored"))
    records += [dict(tool="Python", role="interpreter", version=platform.python_version(), source="https://www.python.org/", source_commit="", licence="PSF-2.0", licence_url_or_source="https://docs.python.org/3/license.html", licence_checked_date="2026-10-03", execution="local", redistributable="yes", notes="project venv"),
                dict(tool="Bash", role="orchestration", version=os.environ["PIPELINE_BASH_VERSION"], source="https://www.gnu.org/software/bash/", source_commit="", licence="GPL-3.0-or-later", licence_url_or_source="https://www.gnu.org/software/bash/", licence_checked_date="2026-10-03", execution="local", redistributable="yes", notes="separately installed executable")]
    table(ROOT / "results/local_software.tsv", records)
    atomic_text(ROOT / "results/system.json", json.dumps(dict(os=platform.platform(), cpu=platform.processor(),
         physical_ram_bytes=psutil.virtual_memory().total, logical_cpus=psutil.cpu_count(), python=sys.executable,
         timestamp=now(), packages={r["tool"]: r["version"] for r in records}), indent=2))
    print("Recorded installed versions and dependency licence metadata")


if __name__ == "__main__":
    main()
