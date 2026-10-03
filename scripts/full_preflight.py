import argparse
import json
from common import ROOT, rows, failure
import importlib


def main():
    p=argparse.ArgumentParser(description="Full mode requires every approved model/target/split artifact")
    p.parse_args()
    problems=["Neural score/recovery aggregation is not yet validated on real returned artifacts. Stage 12 currently merges local methods only; finish and validate the neural evaluator before enabling full mode. See HANDOVER.md."]
    for row in rows(ROOT/"results/gpu_jobs.tsv"):
        name=f'{row["method"]}_{row["target"]}_{row["strategy"]}'
        job=json.loads((ROOT/row["bundle"]/"job.json").read_text())
        directory=ROOT/f"data/remote/{name}"
        if job["licence_status"]!="approved":
            problems.append(f"{name}: licence blocker, see remote/aptadiff/BLOCKER.md")
        elif not (directory/"run_manifest.json").exists():
            problems.append(f"{name}: missing returned artifacts; bash remote/raptgen/run_remote.sh {row['bundle']} output/{name}")
        else:
            importlib.import_module("11_import_generations").validate_artifact(directory,job)
    if problems:
        err=RuntimeError("Full mode refused:\n"+"\n".join(problems))
        failure("full_preflight","GPU artifacts",err,"prepare hosted jobs; resolve licence; import verified outputs")
        raise err


if __name__=="__main__":
    main()
