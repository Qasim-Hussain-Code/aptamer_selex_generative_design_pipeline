"""Read metadata and text only; no external source code is redistributed."""
import argparse
import json
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


if __name__ == "__main__":
    main()
