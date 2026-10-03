import json
from common import ROOT, rows, guard, table


def main():
    guard(800_000_000)
    records = []
    for line in (ROOT / "config/requirements.txt").read_text().splitlines():
        name, version = line.split("==")
        path = ROOT / f"data/external/pypi/{name}.json"
        if not path.exists():
            raise RuntimeError(f"Inspect official PyPI licence/releases first: {name}")
        data = json.loads(path.read_text())
        if version not in data["releases"]:
            raise ValueError(f"Pinned version not verified available: {line}")
        info = data["info"]
        licence = info.get("license_expression") or info.get("license")
        if not licence:
            raise ValueError(f"No package licence metadata: {name}")
        records.append(dict(tool=name, pinned_version=version, current_release=info["version"],
              metadata_sha256=__import__("common").sha256(path), licence=licence,
              install_source=f"https://pypi.org/project/{name}/{version}/", checked_date="2026-10-03",
              rationale="Pinned available versions used in this run; current release inspected; no GPU packages"))
    table(ROOT / "results/package_preflight.tsv", records)


if __name__ == "__main__":
    main()
