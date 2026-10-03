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
        licence = {"numpy":"BSD-3-Clause plus bundled dependency terms", "scipy":"BSD-3-Clause plus bundled dependency terms", "matplotlib":"PSF-derived Matplotlib licence", "ViennaRNA":"ViennaRNA custom licence"}.get(name, licence.splitlines()[0])
        records.append(dict(tool=name, role="local core", version=dist.version,
             source=f"https://pypi.org/project/{name}/{dist.version}/", source_commit="", licence=licence,
             licence_url_or_source="installed distribution METADATA and licence files", licence_checked_date="2026-10-03",
             execution="local", redistributable="yes" if licence != "unclear" else "unclear", notes="Installed separately; no dependency source vendored"))
    records += [dict(tool="Python", role="interpreter", version=platform.python_version(), source="https://www.python.org/", source_commit="", licence="PSF-2.0", licence_url_or_source="https://docs.python.org/3/license.html", licence_checked_date="2026-10-03", execution="local", redistributable="yes", notes="project venv"),
                dict(tool="Bash", role="orchestration", version=os.environ["PIPELINE_BASH_VERSION"], source="https://www.gnu.org/software/bash/", source_commit="", licence="GPL-3.0-or-later", licence_url_or_source="https://www.gnu.org/software/bash/", licence_checked_date="2026-10-03", execution="local", redistributable="yes", notes="separately installed executable")]
    table(ROOT / "results/local_software.tsv", records)
    freeze = subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True)
    atomic_text(ROOT / "config/requirements-lock.txt",freeze)
    transitives=[]
    for dist in sorted(im.distributions(),key=lambda d:d.metadata["Name"].lower()):
        name=dist.metadata["Name"]
        lic=dist.metadata.get("License-Expression") or dist.metadata.get("License") or ";".join(s for s in dist.metadata.get_all("Classifier",[]) if "License" in s) or "unclear"
        lic={"numpy":"BSD-3-Clause plus bundled dependency terms", "scipy":"BSD-3-Clause plus bundled dependency terms", "matplotlib":"PSF-derived Matplotlib licence", "ViennaRNA":"ViennaRNA custom licence", "cycler":"BSD-3-Clause", "kiwisolver":"BSD-3-Clause", "python-dateutil":"Apache-2.0 and BSD-3-Clause contributions"}.get(name,lic.splitlines()[0][:180])
        evidence=[f for f in dist.files or [] if "licen" in str(f).lower() and str(f).endswith(("LICENSE","LICENSE.txt","LICENSE.rst"))]
        transitives.append(dict(tool=name,version=dist.version,licence=lic,source=f"https://pypi.org/project/{name}/{dist.version}/",installed_separately="yes",checked_date="2026-10-03",licence_evidence=";".join(str(f) for f in evidence)))
    table(ROOT/"results/transitive_licences.tsv",transitives)
    atomic_text(ROOT / "results/system.json", json.dumps(dict(os=platform.platform(), cpu=platform.processor(),
         physical_ram_bytes=psutil.virtual_memory().total, logical_cpus=psutil.cpu_count(), python=sys.executable,
         timestamp=now(), packages={r["tool"]: r["version"] for r in records}), indent=2))
    print("Recorded installed versions and dependency licence metadata")


if __name__ == "__main__":
    main()
