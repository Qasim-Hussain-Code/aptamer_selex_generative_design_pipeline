"""Acquire small reference inputs independently so a failed service cannot block others."""
import argparse
import json
import urllib.parse
import zipfile
from common import ROOT, download, failure


def fetch(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    try:
        return download(url, dest)
    except Exception as e:
        failure("references", url, e, "retain independent inputs; retry exact URL")
        print(type(e).__name__, url, str(e)[:200], flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--part", choices=["metadata", "upstream", "packages", "supplement"], required=True)
    a = p.parse_args()
    if a.part == "metadata":
        fields = "run_accession,experiment_accession,sample_accession,study_accession,study_title,experiment_title,sample_title,fastq_ftp,fastq_md5,fastq_bytes,read_count"
        for acc in ["DRA009383", "DRA009384"]:
            query = urllib.parse.urlencode(dict(result="read_run", query=f'submission_accession="{acc}"', fields=fields, format="tsv"))
            fetch("https://www.ebi.ac.uk/ena/portal/api/search?" + query, ROOT / f"data/external/{acc}.ena.tsv")
            # The direct DDBJ route failed during source inspection. Preserve that failure
            # record and use the official INSDC mirror rather than retrying a guessed path.
    elif a.part == "supplement":
        url = "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/48/14/10.1093_nar_gkaa484/2/gkaa484_supplemental_file.pdf?Expires=2147483647&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA&Signature=oftgbjSTc2kY8tnC9FlAXtizeAEfhTN1Rh-rBDyk~-Gz~u3VsEBszccmYbZDXh7MudW4vuml3eP~RjGVSoSvcyvB~e2eSt7HrKmbIw7nPhso4SW8YNjl8QKeor5z~QGwrGC~IfjiK3E9Iq0pnLLnrDo8b55kT0xHbOF7kPWlvnkqlRxvDxENWu5Kn9KZe0~olRbLr3Joneom~ZHHRwN1u3vcbhZV2F9VINQY~qsx2fMp6TKdxx5nW8n9e0JUrbmCNbTIBH7Us7EX~Bke7aacBm~q4GHhqT0O-VpJJ~g9pCeLD~00YTpBG2Yb0ptiHRPlWHmf4XwJAHPp6ov-WhkutQ__"
        fetch(url, ROOT / "data/external/raptranker_supplement.pdf")
    elif a.part == "upstream":
        for name in ["raptgen", "aptadiff", "instructna"]:
            d = ROOT / f"data/external/{name}"
            info = json.loads((d / "repository.json").read_text())
            files = json.loads((d / "tree.json").read_text())["tree"]
            commit = json.loads((d / "commits.json").read_text())[0]["sha"]
            print(name, commit, info["license"], flush=True)
            for f in files:
                path = f["path"]
                wanted = path in ["LICENSE", "README.md", "Dockerfile", "Pipfile", "requirements.txt", "scripts/real.py", "scripts/decode.py", "scripts/gmm.py"]
                wanted |= name == "aptadiff" and path.startswith("data/raw_data/")
                if f["type"] == "blob" and wanted:
                    fetch(f'https://raw.githubusercontent.com/{info["full_name"]}/{commit}/{path}', d / path)
    elif a.part == "packages":
        for pkg in ["ViennaRNA", "numpy", "scipy", "scikit-learn", "matplotlib", "PyYAML", "pytest", "psutil", "pypdf", "rapidfuzz", "shellcheck_py"]:
            dest = ROOT / f"data/external/pypi/{pkg}.json"
            if fetch(f"https://pypi.org/pypi/{pkg}/json", dest):
                info = json.loads(dest.read_text())["info"]
                print(pkg, info["version"], info.get("license_expression"), str(info.get("license"))[:80], flush=True)


if __name__ == "__main__":
    main()
