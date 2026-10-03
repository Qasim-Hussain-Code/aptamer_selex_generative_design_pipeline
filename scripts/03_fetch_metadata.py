"""Resolve rounds from experiment library names, never from accession order."""
import argparse
import re
import xml.etree.ElementTree as ET
from common import ROOT, rows, table, download, failure


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", choices=["all", "tg2", "integrin"], default="all")
    p.add_argument("--offline", action="store_true")
    a = p.parse_args()
    out = []
    for d in rows(ROOT / "config/datasets.tsv"):
        if a.target != "all" and a.target != d["target"]:
            continue
        acc = d["accession"]
        source = ROOT / f"data/external/{acc}.ena.tsv"
        if not source.exists():
            raise RuntimeError(f"Missing ENA metadata: bash scripts/acquire_reference_inputs.sh --part metadata ({acc})")
        for r in rows(source):
            exp = r["experiment_accession"]
            xml = ROOT / f"data/external/experiments/{exp}.xml"
            if not xml.exists() and not a.offline:
                download(f"https://www.ebi.ac.uk/ena/browser/api/xml/{exp}", xml, exp)
            if not xml.exists():
                raise RuntimeError(f"Experiment record not available: {exp}")
            tree = ET.parse(xml)
            title = " ".join(x.text or "" for x in tree.iter() if x.tag in ["TITLE", "LIBRARY_NAME"])
            aliases = " ".join(x.attrib.get("alias", "") for x in tree.iter())
            text = title + " " + aliases
            matches = re.findall(r"(?:round[_ -]?(\d+)|(\d+)[Rr])", text, re.I)
            rounds = {int(x or y) for x, y in matches}
            if len(rounds) != 1:
                raise ValueError(f"Cannot unambiguously assign round to {exp}: {text}")
            rnd = rounds.pop()
            expected = [int(x) for x in d["rounds"].split(",")]
            if rnd not in expected:
                raise ValueError(f"Unexpected round {rnd}: {exp}")
            for url, md5, size in zip(r["fastq_ftp"].split(";"), r["fastq_md5"].split(";"), r["fastq_bytes"].split(";"), strict=True):
                out.append(dict(target=d["target"], accession=acc, run=r["run_accession"], experiment=exp,
                    study=r["study_accession"], study_title=r["study_title"], round=rnd, bytes=int(size),
                    md5=md5, url="https://" + url, expected_reads=r["read_count"],
                    round_evidence=text, metadata_source=f"https://www.ebi.ac.uk/ena/browser/api/xml/{exp}"))
    out.sort(key=lambda x: (x["target"], x["round"], x["run"]))
    for target in {r["target"] for r in out}:
        rounds_found = [r["round"] for r in out if r["target"] == target]
        if len(set(rounds_found)) != len(rounds_found):
            raise ValueError("Multiple files per round require explicit aggregation; refusing")
    table(ROOT / "results/selex_manifest.tsv", out)
    print(f"Validated metadata for {len(out)} FASTQ files; no sequence download yet")


if __name__ == "__main__":
    main()
