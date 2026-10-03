import importlib
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from common import table,sha256
from configure import effective_budget


@pytest.mark.parametrize("requested,free,expected",[(13_000_000_000,14_100_000_000,13_000_000_000),(20_000_000_000,30_000_000_000,13_000_000_000),(13_000_000_000,5_000_000_000,4_000_000_000)])
def test_decimal_disk_budget(requested,free,expected):
    assert effective_budget(requested,free)==expected


def test_reserve_refusal():
    with pytest.raises(ValueError):
        effective_budget(13_000_000_000,900_000_000)


def test_projected_peak_refusal(tmp_path,monkeypatch):
    import common
    monkeypatch.setattr(common,"disk_bytes",lambda root:100)
    monkeypatch.setattr(common.shutil,"disk_usage",lambda root:SimpleNamespace(free=10000))
    monkeypatch.setenv("DISK_BUDGET_BYTES","1000")
    monkeypatch.setenv("RESERVE_BYTES","100")
    assert common.guard(500,tmp_path)==650
    with pytest.raises(RuntimeError,match="shortfall"):
        common.guard(900,tmp_path)


def test_checkpoint_hash(tmp_path):
    module=importlib.import_module("11_import_generations")
    path=tmp_path/"checkpoint"
    path.write_bytes(b"known-checkpoint-fixture")
    module.verify_hash(path,sha256(path))
    with pytest.raises(ValueError):
        module.verify_hash(path,"0"*64)


@pytest.fixture
def artifact(tmp_path):
    job=dict(method="raptgen",target="tg2",strategy="family",source_repository="https://github.com/hmdlab/raptgen",source_commit="c4986ca9fa439b9389916c05829da4ff9c30d6f3",configuration_hash="a"*64,input_hashes={"train.tsv":"b"*64},context="unsupervised_or_selex_only_generation",licence_status="approved",generation_seeds=[11],training_seed=1729,generation_budget=2,minimum_length=4,maximum_length=4,forbidden_primers=["CCCC"])
    table(tmp_path/"generations_11.tsv.gz",[dict(sequence="ACGT"),dict(sequence="TGCA")])
    table(tmp_path/"heldout_scores.tsv",[dict(sequence="ACGT",score=-1,score_type="posterior_mean_conditional_PHMM_log_probability_per_nt")])
    manifest={k:job[k] for k in ["method","target","strategy","source_repository","source_commit","configuration_hash","input_hashes","context","generation_seeds","training_seed"]}
    manifest.update(status="complete",assay_values_used=False,package_versions=["torch==1.5.0"],cuda_version="10.2",gpu_model="synthetic_fixture",gpu_memory_bytes=1,peak_allocated_gpu_memory_bytes=1,start_time="2026-10-03T00:00:00Z",end_time="2026-10-03T00:00:01Z",elapsed_seconds=1,checkpoint_hash="c"*64,generation_statistics=[dict(seed=11,raw=2,valid=2,unique=2)],output_hashes={n:sha256(tmp_path/n) for n in ["generations_11.tsv.gz","heldout_scores.tsv"]})
    (tmp_path/"run_manifest.json").write_text(json.dumps(manifest))
    return tmp_path,job,manifest


def test_imported_artifact(artifact):
    path,job,manifest=artifact
    assert importlib.import_module("11_import_generations").validate_artifact(path,job)["status"]=="complete"


@pytest.mark.parametrize("field,value",[("source_commit","wrong"),("configuration_hash","wrong"),("assay_values_used",True),("generation_seeds",[99]),("checkpoint_hash","not_a_hash")])
def test_manifest_tampering_refused(artifact,field,value):
    path,job,manifest=artifact
    manifest[field]=value
    (path/"run_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        importlib.import_module("11_import_generations").validate_artifact(path,job)


def test_output_corruption_refused(artifact):
    path,job,manifest=artifact
    (path/"heldout_scores.tsv").write_text("tampered")
    with pytest.raises(ValueError):
        importlib.import_module("11_import_generations").validate_artifact(path,job)


def test_unclear_licence_refused(artifact):
    path,job,manifest=artifact
    job["licence_status"]="blocked"
    with pytest.raises(ValueError,match="licence"):
        importlib.import_module("11_import_generations").validate_artifact(path,job)


def test_primary_pdf_parser_missing_rows_refused():
    with pytest.raises(ValueError):
        importlib.import_module("06_build_ground_truth").parse_text("Data1-1 ACGT 161007 340.2 1")


@pytest.mark.parametrize("changed",["input","output","code","configuration"])
def test_restart_rejects_changed_content(tmp_path,changed,monkeypatch):
    import common
    monkeypatch.delenv("MASTER_SEED",raising=False)
    (tmp_path/"config").mkdir()
    (tmp_path/"scripts").mkdir()
    config=tmp_path/"config/benchmark.yml"
    config.write_text("master_seed: 1729\n")
    code=tmp_path/"scripts/common.py"
    code.write_text("first code version\n")
    source=tmp_path/"input.tsv";source.write_text("input\n")
    output=tmp_path/"output.tsv";output.write_text("output\n")
    common.stamp("fixture",[source],[output],root=tmp_path)
    assert common.completed("fixture",[source],[output],root=tmp_path)
    paths={"input":source,"output":output,"code":code,"configuration":config}
    paths[changed].write_text("master_seed: 99\n" if changed=="configuration" else "changed content\n")
    assert not common.completed("fixture",[source],[output],root=tmp_path)
