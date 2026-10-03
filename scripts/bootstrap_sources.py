import argparse
import json
from common import ROOT, download


def main():
    p=argparse.ArgumentParser(description="Fetch upstream inventory and primary XML only when missing")
    p.parse_args()
    for name,repo in [("raptgen","hmdlab/raptgen"),("aptadiff","wz-create/AptaDiff"),("instructna","zhimingzhang275/InstructNA")]:
        for suffix,filename in [("","repository.json"),("/commits?per_page=1","commits.json"),("/git/trees/HEAD?recursive=1","tree.json")]:
            dest=ROOT/f"data/external/{name}/{filename}"
            if not dest.exists():
                download("https://api.github.com/repos/"+repo+suffix,dest)
    for name,pmc in [("raptranker","PMC7641312"),("aptadiff","PMC11491854")]:
        dest=ROOT/f"data/external/{name}.xml"
        if not dest.exists():
            download(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML",dest)


if __name__=="__main__":
    main()
