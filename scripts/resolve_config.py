import argparse
import hashlib
import json
import yaml
from common import ROOT, settings, atomic_text, table, now


def main():
    p=argparse.ArgumentParser(description="Record immutable resolved configuration; archive changes")
    p.parse_args()
    cfg=settings()
    text=yaml.safe_dump(cfg,sort_keys=True)
    out=ROOT/"results/resolved_config.yml"
    if out.exists() and out.read_text()!=text:
        old_hash=hashlib.sha256(out.read_bytes()).hexdigest()
        atomic_text(ROOT/f"logs/config_history/{old_hash}.yml",out.read_text())
        print("Previous resolved configuration preserved in logs/config_history")
    atomic_text(out,text)
    table(ROOT/"results/configuration_hash.tsv",[dict(configuration_hash=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest(),resolved_file_sha256=hashlib.sha256(text.encode()).hexdigest(), master_seed=cfg["master_seed"], recorded_at=now())])


if __name__=="__main__":
    main()
