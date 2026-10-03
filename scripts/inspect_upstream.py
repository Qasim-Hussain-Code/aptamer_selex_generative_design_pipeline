"""Read metadata and text only; no external source code is redistributed."""
import argparse
import json
import urllib.parse
import xml.etree.ElementTree as ET
from common import ROOT, download, atomic_text


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    for name, repo in [("raptgen", "hmdlab/raptgen"), ("aptadiff", "wz-create/AptaDiff"), ("instructna", "zhimingzhang275/InstructNA")]:
        base = f"https://api.github.com/repos/{repo}"
        for suffix, filename in [("", "repository.json"), ("/commits?per_page=1", "commits.json"), ("/git/trees/HEAD?recursive=1", "tree.json")]:
            download(base + suffix, ROOT / f"data/external/{name}/{filename}")
    for accession in ["DRA009383", "DRA009384"]:
        download(f"https://ddbj.nig.ac.jp/resource/sra-submission/{accession}.xml", ROOT / f"data/external/{accession}.xml", accession)
    download("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7641312/fullTextXML", ROOT / "data/external/raptranker.xml")
    download("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11491854/fullTextXML", ROOT / "data/external/aptadiff.xml")
    download("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7641312/supplementaryFiles", ROOT / "data/external/raptranker_supplements.zip")
    fields = "run_accession,experiment_accession,sample_accession,study_accession,study_title,experiment_title,sample_title,fastq_ftp,fastq_md5,fastq_bytes,read_count"
    for accession in ["DRA009383", "DRA009384"]:
        query = urllib.parse.urlencode(dict(result="read_run", query=f'submission_accession="{accession}"', fields=fields, format="tsv"))
        download("https://www.ebi.ac.uk/ena/portal/api/search?" + query, ROOT / f"data/external/{accession}.ena.tsv", accession)
    for name in ["raptgen", "aptadiff", "instructna"]:
        d = ROOT / f"data/external/{name}"
        info = json.loads((d / "repository.json").read_text())
        tree = json.loads((d / "tree.json").read_text())["tree"]
        commit = json.loads((d / "commits.json").read_text())[0]["sha"]
        print(name, commit, "license", info["license"])
        files = [f for f in tree if f["type"] == "blob"]
        print("FILES", [(f["path"], f.get("size")) for f in files if name != "instructna" or "licen" in f["path"].lower()])
        wanted = [f for f in files if f["path"] in ["LICENSE", "README.md", "Dockerfile", "Pipfile", "requirements.txt", "scripts/real.py", "scripts/gmm.py", "scripts/sample.py"]]
        for f in wanted:
            download(f'https://raw.githubusercontent.com/{info["full_name"]}/{commit}/{f["path"]}', d / f["path"])
    for pkg in ["ViennaRNA", "numpy", "scipy", "scikit-learn", "matplotlib", "PyYAML", "pytest", "psutil", "pypdf", "rapidfuzz"]:
        dest = ROOT / f"data/external/pypi/{pkg}.json"
        download(f"https://pypi.org/pypi/{pkg}/json", dest)
        info = json.loads(dest.read_text())["info"]
        print("PYPI", pkg, info["version"], info.get("license_expression"), str(info.get("license"))[:100], info["project_urls"])


if __name__ == "__main__":
    main()
