"""Write the scientific narrative only from completed result tables and logs."""
import argparse
from collections import defaultdict
import html
import json
import re
from pathlib import Path
from datetime import datetime, timedelta, timezone
from common import ROOT, rows, table, atomic_text, settings, sha256


def markdown_table(headers, data):
    return "| " + " | ".join(headers) + " |\n| " + " | ".join(["---"]*len(headers)) + " |\n" + "\n".join("| "+" | ".join(map(str,r))+" |" for r in data)


def render_template(template, values):
    missing=set(re.findall(r"\{\{([a-z0-9_]+)\}\}",template))-set(values)
    if missing:
        raise ValueError("Missing README values: "+", ".join(sorted(missing)))
    for key,value in values.items():
        template=template.replace("{{"+key+"}}",str(value))
    if "{{" in template or "}}" in template:
        raise ValueError("Unresolved README template placeholder")
    return template


def report_html(markdown):
    def href(url):
        return url if url.startswith(("https://","http://")) else "../"+url
    def inline(value):
        value=html.escape(value)
        value=re.sub(r'!\[([^]]*)\]\(([^)]+)\)',lambda m:f'<img src="{href(m[2])}" alt="{m[1]}">',value)
        value=re.sub(r'(?<!!)\[([^]]*)\]\(([^)]+)\)',lambda m:f'<a href="{href(m[2])}">{m[1]}</a>',value)
        value=re.sub(r'`([^`]+)`',r'<code>\1</code>',value)
        return re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',value)
    output=[];lines=markdown.splitlines();i=0
    while i<len(lines):
        line=lines[i]
        if not line.strip():i+=1;continue
        if line.startswith("```"):
            block=[];i+=1
            while i<len(lines) and not lines[i].startswith("```"):
                block.append(lines[i]);i+=1
            output.append("<pre><code>"+html.escape("\n".join(block))+"</code></pre>");i+=1;continue
        if line.startswith("|"):
            block=[]
            while i<len(lines) and lines[i].startswith("|"):
                block.append([c.strip() for c in lines[i].strip("|").split("|")]);i+=1
            output.append('<div class="table"><table><thead><tr>'+''.join('<th>'+inline(c)+'</th>' for c in block[0])+"</tr></thead><tbody>"+''.join('<tr>'+''.join('<td>'+inline(c)+'</td>' for c in row)+'</tr>' for row in block[2:])+"</tbody></table></div>");continue
        heading=re.match(r"^(#{1,6}) (.+)",line)
        if heading:
            level=len(heading[1]);output.append(f"<h{level}>"+inline(heading[2])+f"</h{level}>");i+=1;continue
        if line.startswith("- "):
            items=[]
            while i<len(lines) and lines[i].startswith("- "):
                items.append("<li>"+inline(lines[i][2:])+"</li>");i+=1
            output.append("<ul>"+"".join(items)+"</ul>");continue
        paragraph=[line];i+=1
        while i<len(lines) and lines[i].strip() and not lines[i].startswith(("#","|","```","- ")):
            paragraph.append(lines[i]);i+=1
        output.append("<p>"+inline(" ".join(paragraph))+"</p>")
    return '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Executed local analysis</title><style>body{max-width:1100px;margin:40px auto;padding:0 20px;font:16px/1.6 system-ui;color:#182c3b}img{max-width:100%}pre{overflow:auto;background:#edf2f5;padding:16px}code{font-size:13px}.table{overflow:auto}table{border-collapse:collapse;width:100%}th,td{text-align:left;border-bottom:1px solid #ccd8df;padding:8px}h2{margin-top:40px}a{color:#286181}</style>'+"\n".join(output)+"</html>"


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    cfg=settings()
    metrics=list(rows(ROOT/"results/evaluation_metrics.tsv"))
    intervals=list(rows(ROOT/"results/bootstrap_intervals.tsv"))
    paired=list(rows(ROOT/"results/paired_differences.tsv"))
    predictions=list(rows(ROOT/"results/held_out_predictions.tsv"))
    qc=list(rows(ROOT/"results/preprocessing_counts.tsv"))
    truth_audit={r["target"]:r for r in rows(ROOT/"results/ground_truth_audit.tsv")}
    leakage=list(rows(ROOT/"results/leakage_audit.tsv"))
    generation=list(rows(ROOT/"results/generation_statistics.tsv"))
    resources=list(rows(ROOT/"logs/resource_usage.tsv"))
    software=list(rows(ROOT/"results/software_manifest.tsv"))
    primary=cfg["family_threshold"]
    def metric(t,m,strategy="family",name="average_precision"):
        return float(next(r["value"] for r in metrics if r["target"]==t and r["method"]==m and r["strategy"]==strategy and float(r["threshold"])==primary and r["metric"]==name))
    def interval(t,m):
        r=next(r for r in intervals if r["target"]==t and r["method"]==m and r["strategy"]=="family" and float(r["threshold"])==primary)
        return f'{float(r["lower_95"]):.3f} to {float(r["upper_95"]):.3f}'
    def delta(t):
        return next(r for r in paired if r["target"]==t and r["method"]=="kmer_structure_ridge" and r["strategy"]=="family" and float(r["threshold"])==primary)
    peak=max(int(r["peak_observed_disk_bytes"]) for r in resources)
    rss=max(int(r["peak_rss_bytes"]) for r in resources)
    minfree=min(int(r["minimum_free_bytes"]) for r in resources)
    total_candidates=sum(int(r["unique_generations"]) for r in generation)
    exact_recovery=sum(int(r["exact_test_recovery"]) for r in generation)
    family_recovery=sum(int(r["positive_test_family_recovery"]) for r in generation)
    summary=[]
    for target in ["tg2","integrin"]:
        data=[r for r in qc if r["target"]==target]
        aud=next(r for r in leakage if r["target"]==target and r["strategy"]=="family" and float(r["threshold"])==primary)
        test=[r for r in predictions if r["target"]==target and r["method"]=="kmer_ridge" and r["strategy"]=="family" and float(r["threshold"])==primary]
        summary.append(dict(target=target,raw_reads=sum(int(r["raw_reads"]) for r in data),retained_reads=sum(int(r["retained_reads"]) for r in data),
             discarded_reads=sum(int(r["discarded_reads"]) for r in data),rounds=len(data),test_sequences=len(test),test_positives=sum(int(r["experimental_label"]) for r in test),test_negatives=sum(1-int(r["experimental_label"]) for r in test),
             test_families=aud["test_families"],training_sequences=aud["training_sequences"],selected_pool=aud["selected_pool"],pool_cap_exclusions=aud["excluded_from_pool"],
             family_exclusions=aud["excluded_training_sequences"],ridge_ap=metric(target,"kmer_ridge"),observed_enrichment_ap=metric(target,"enrichment_observed"),
             assays_retained=truth_audit[target]["retained"],assays_excluded=truth_audit[target]["excluded"],
             peak_project_bytes=peak,peak_measured_rss_bytes=rss,minimum_free_bytes=minfree,resource_snapshot_rows=len(resources)))
    table(ROOT/"results/summary.tsv",summary)
    table(ROOT/"results/dataset_design.tsv",list(rows(ROOT/"config/datasets.tsv")))
    run_summary=dict(run="executed_core",unique_markov_candidates=total_candidates,generation_runs=len(generation),assay_records_retained=sum(int(r["retained"]) for r in truth_audit.values()),exact_test_recovery=exact_recovery,positive_test_family_recovery=family_recovery,
                     primary_family_threshold=primary,bootstrap_replicates=cfg["bootstrap_replicates"],generation_budget=cfg["generation_budget"],
                     generation_seeds=",".join(map(str,cfg["generation_seeds"])),folding_temperature_c=cfg["folding_temperature_c"],resource_snapshot_rows=len(resources))
    table(ROOT/"results/run_summary.tsv",[run_summary])
    sources=[]
    for r in summary:
        for k,v in r.items():
            if k!="target": sources.append(dict(claim_id=f'{r["target"]}_{k}',value=v,source_file="results/summary.tsv",source_row=f'target={r["target"]}',source_column=k,scope="measured_or_manifest_count"))
    for r in metrics:
        if float(r["threshold"])==primary:
            sources.append(dict(claim_id=f'{r["target"]}_{r["strategy"]}_{r["method"]}_{r["metric"]}',value=r["value"],source_file="results/evaluation_metrics.tsv",source_row=f'target={r["target"]};strategy={r["strategy"]};method={r["method"]};metric={r["metric"]};threshold={primary}',source_column="value",scope="measured"))
    for r in paired:
        if float(r["threshold"])==primary:
            for col in ["difference","lower_95","upper_95"]:
                sources.append(dict(claim_id=f'{r["target"]}_{r["strategy"]}_{r["method"]}_{r["comparison"]}_{col}',value=r[col],source_file="results/paired_differences.tsv",source_row=f'target={r["target"]};strategy={r["strategy"]};method={r["method"]};comparison={r["comparison"]};threshold={primary}',source_column=col,scope="measured"))
    for k,v in run_summary.items():
        if k!="run":sources.append(dict(claim_id=k,value=v,source_file="results/run_summary.tsv",source_row="run=executed_core",source_column=k,scope="aggregate_or_configured"))
    for r in intervals:
        if float(r["threshold"])==primary:
            for col in ["lower_95","upper_95"]:
                sources.append(dict(claim_id=f'{r["target"]}_{r["strategy"]}_{r["method"]}_ap_{col}',value=r[col],source_file="results/bootstrap_intervals.tsv",source_row=f'target={r["target"]};strategy={r["strategy"]};method={r["method"]};threshold={primary}',source_column=col,scope="measured_family_bootstrap"))
    latest_resources={r["stage"]:r for r in resources}
    for stage,r in latest_resources.items():
        for col in ["elapsed_seconds","peak_rss_bytes","peak_observed_disk_bytes"]:
            sources.append(dict(claim_id=f'resource_{stage}_{col}',value=r[col],source_file="logs/resource_usage.tsv",source_row=f'stage={stage};start_timestamp={r["start_timestamp"]}',source_column=col,scope="raw_units; prose_MB_and_GB_use_decimal_divisors"))
    table(ROOT/"results/readme_traceability.tsv",sources)
    definitions=[
        ("relative_frequency","count/retained_round_depth","fraction","higher","selection readout, not experimental affinity"),
        ("log2_enrichment","log2((c_late+a)/(N_late+aV) / ((c_previous+a)/(N_previous+aV)))","dimensionless","higher","amplification and bottleneck biases; a from config"),
        ("predicted_log2_enrichment","regularized ridge prediction using training SELEX enrichment","dimensionless","higher","not SPR or KD prediction"),
        ("decision_function","regularized logistic decision function fitted on training published labels","dimensionless","higher","uncalibrated ranking score"),
        ("per_nt_log_probability","first-order Markov natural log sequence probability / sequence length","nats/nt","higher","training distribution likelihood, not activity"),
        ("posterior_mean_conditional_PHMM_log_probability_per_nt","RaptGen conditional reconstruction PHMM log probability / length","nats/nt","higher","registered but no neural scores generated"),
        ("AptaDiff","unregistered; method excluded","not_available","not_available","no diffusion objective is interpreted as activity"),
        ("SPR_response","source end-of-injection SPR response","RU","higher","experimental response, not KD"),
        ("published_label","original Positive/Negative (1/0)","binary","1_positive","published designation preserved"),
        ("KD_strength","-log10(KD)","log_units","higher","direction tested; no KD used in current evaluation"),
        ("mfe_kcal_mol","ViennaRNA intramolecular MFE","kcal/mol","no_binding_direction","canonical-RNA secondary-structure proxy"),
        ("base_pair_count","opening parentheses in MFE dot-bracket","pairs","no_binding_direction","intramolecular structure"),
        ("paired_fraction","2*base_pair_count/construct_length","fraction","no_binding_direction","full assay construct excluding polyA tether"),
        ("ensemble_diversity","partition-function mean base-pair distance","pairs","no_binding_direction","canonical RNA ensemble proxy"),
        ("average_precision","sum precision at distinct score thresholds * recall increment","fraction","higher","primary; compare prevalence and test class count"),
        ("auroc","probability a positive is ranked above a negative with half credit to ties","fraction","higher","secondary; undefined with one class"),
        ("spearman","rank correlation of model score and source RU","correlation","higher","undefined for constant scores"),
        ("kendall","Kendall tau of model score and source RU","correlation","higher","secondary; ties included"),
        ("precision_top_k","positive fraction among top k with fractional tied boundary","fraction","higher","k fixed before evaluation"),
        ("prevalence","test positives / test sequences","fraction","not_a_model_score","constant-score AP reference"),
        ("maximum_identity","max training (1 - global edit distance / max lengths)","fraction","no_binding_direction","not an alignment matches percentage"),
        ("edit_distance","minimum global Levenshtein edits to training","edits","no_binding_direction","insertions and deletions included"),
        ("family_recovery","test family with identity edge from at least one generated candidate","families","higher_coverage","no untested candidate is established as a binder"),
        ("gc_fraction","G+C nucleotide count / total nucleotides","fraction","no_binding_direction","generation composition"),
        ("duplicate_fraction","1 - unique/retained reads","fraction","not_a_binding_score","within-round count redundancy"),
        ("paired_difference","metric_left - metric_right on identical test sequences","metric_units","context_specific","family bootstrap; excludes training uncertainty"),
        ("seed_standard_deviation","sample standard deviation across all fixed generation seeds","quantity_units","not_a_binding_score","does not estimate training-seed variance")]
    table(ROOT/"results/score_dictionary.tsv",[dict(score=n,definition=d,units=u,direction=direction,limitations=limit) for n,d,u,direction,limit in definitions])
    st={r["target"]:r for r in summary}
    names={"tg2":"TG2", "integrin":"Integrin alpha V beta 3"}
    design={r["target"]:r for r in rows(ROOT/"config/datasets.tsv")}
    system=json.loads((ROOT/"results/system.json").read_text())
    downloads=list(rows(ROOT/"results/download_manifest.tsv"))
    local_dates=sorted({datetime.fromisoformat(r["retrieval_timestamp"]).astimezone(timezone(timedelta(hours=8))).date().isoformat() for r in downloads})
    retrieval_date=(local_dates[0] if len(local_dates)==1 else local_dates[0]+" to "+local_dates[-1])+" (Asia/Taipei)"
    provenance=dict(scope="readme_sources",retrieval_date_taipei=retrieval_date,download_manifest_sha256=sha256(ROOT/"results/download_manifest.tsv"),
                    python_version=system["packages"]["Python"],bash_version=system["packages"]["Bash"],viennarna_version=system["packages"]["ViennaRNA"])
    table(ROOT/"results/readme_provenance.tsv",[provenance])
    primary_ci={t:next(r for r in intervals if r["target"]==t and r["strategy"]=="family" and r["method"]=="kmer_ridge" and float(r["threshold"])==primary) for t in st}
    split_delta={t:next(r for r in paired if r["target"]==t and r["strategy"]=="random" and r["method"]=="kmer_ridge" and r["comparison"]=="random_minus_family" and float(r["threshold"])==primary) for t in st}
    late_qc={t:max((r for r in qc if r["target"]==t),key=lambda r:int(r["round"])) for t in st}
    method_labels={"frequency":"Training frequency lookup", "enrichment":"Training enrichment lookup", "kmer_ridge":"SELEX enrichment ridge", "markov":"Markov log probability per nt", "kmer_logistic":"Assay-label logistic", "kmer_structure_ridge":"SELEX ridge + structure", "kmer_structure_logistic":"Assay logistic + structure", "frequency_observed":"Observed frequency", "enrichment_observed":"Observed enrichment"}
    ap_table=markdown_table(["Method","TG2 AP (95% interval)","Integrin AP (95% interval)"],[[label,f'{metric("tg2",m):.3f} ({interval("tg2",m)})',f'{metric("integrin",m):.3f} ({interval("integrin",m)})'] for m,label in method_labels.items()])
    data_table=markdown_table(["Target / accession","Rounds","Raw reads","Retained reads","Rejected reads","Assays retained / excluded"],[[f'{names[t]} / {design[t]["accession"]}',design[t]["rounds"].split(",")[0]+"-"+design[t]["rounds"].split(",")[-1],f'{r["raw_reads"]:,}',f'{r["retained_reads"]:,}',f'{r["discarded_reads"]:,}',f'{r["assays_retained"]} / {r["assays_excluded"]}'] for t,r in st.items()])
    test_table=markdown_table(["Primary panel","Test sequences","Positive / negative","Independent test families"],[[names[t],r["test_sequences"],f'{r["test_positives"]} / {r["test_negatives"]}',r["test_families"]] for t,r in st.items()])
    split_labels={"naive_sequence_random":"Independent sequence partition", "random":"Matched exact-only exclusion", "family":"Family exclusion"}
    split_rows=[]
    for t in st:
        for strategy,label in split_labels.items():
            record=next(r for r in metrics if r["target"]==t and r["strategy"]==strategy and r["method"]=="kmer_ridge" and r["metric"]=="average_precision" and float(r["threshold"])==primary)
            split_rows.append([names[t],label,record["test_sequences"],f'{metric(t,"kmer_ridge",strategy,"prevalence"):.3f}',f'{float(record["value"]):.3f}'])
            sources.append(dict(claim_id=f'{t}_{strategy}_test_count',value=record["test_sequences"],source_file="results/evaluation_metrics.tsv",source_row=f'target={t};strategy={strategy};method=kmer_ridge;metric=average_precision;threshold={primary}',source_column="test_sequences",scope="measured"))
    split_table=markdown_table(["Target","Split","Test sequences","Positive fraction","Ridge AP"],split_rows)
    executed=[]
    for stage,state_name in [("splits","splits"),("baselines","baselines"),("secondary","secondary"),("ablation","ablation"),("candidates","generation"),("evaluation","evaluation")]:
        state_path=ROOT/f"logs/state/{state_name}.json"
        if not state_path.exists():continue
        updated=datetime.fromisoformat(json.loads(state_path.read_text())["time"])
        matches=[r for r in resources if r["stage"]==stage and r["exit_status"]=="0" and datetime.fromisoformat(r["start_timestamp"])<=updated<=datetime.fromisoformat(r["end_timestamp"])]
        if not matches:continue
        r=matches[-1]
        executed.append(dict(stage=stage,start_timestamp=r["start_timestamp"],elapsed_seconds=r["elapsed_seconds"],peak_rss_bytes=r["peak_rss_bytes"],peak_observed_disk_bytes=r["peak_observed_disk_bytes"],selection="successful_execution_containing_scientific_manifest_update"))
    table(ROOT/"results/readme_resource_summary.tsv",executed)
    resource_table=markdown_table(["Scientific stage","Elapsed s","Peak sampled RSS MB","Peak sampled project GB"],[[r["stage"],r["elapsed_seconds"],f'{int(r["peak_rss_bytes"])/1e6:.1f}',f'{int(r["peak_observed_disk_bytes"])/1e9:.3f}'] for r in executed])
    stage_definitions=[
        ("00 configure","configure","","bash scripts/00_configure.sh --threads 1 --ram 16 --disk 13 --gpu-mode none --yes","Host resources -> project.conf; enforce memory, disk and thread limits."),
        ("01 sources","sources","sources","stage sources scripts/01_verify_sources.py --offline","Acquired publications and repository metadata -> provenance and licence records; establish source identity."),
        ("02 install","install","","bash scripts/02_install.sh","Pinned package versions -> CPU environment; keep dependency and chemistry assumptions explicit."),
        ("03 metadata","metadata","","stage metadata scripts/03_fetch_metadata.py","Accessions and experiment XML -> run manifest; resolve rounds from source evidence."),
        ("04/05 reads","fetch","","bash scripts/04_fetch_selex.sh --target all","One FASTQ at a time -> verified compact counts and QC; displayed time is cached validation, with per-round processing listed separately."),
        ("06 assays","ground_truth","ground_truth","stage ground_truth scripts/06_build_ground_truth.py","Supplement S4/S5 -> assay manifest; preserve original labels and SPR units."),
        ("07 splits","splits","splits","stage splits scripts/07_build_splits.py","Late-round pool and assays -> family partitions; exclude test relatives before fitting."),
        ("08 baselines","baselines","baselines","stage baselines scripts/08_baselines.py","Training counts and permitted assays -> fixed sequence fits and scores; retain exposure-matched baselines."),
        ("09 structure","secondary","secondary","stage secondary scripts/09_secondary_structure.py","Complete assay constructs -> canonical-RNA descriptors; approximate intramolecular structure."),
        ("10 GPU jobs","gpu_jobs","","stage gpu_jobs scripts/10_prepare_gpu_jobs.py","Training-only exports -> pinned RaptGen job bundles; isolate held-out assay responses."),
        ("11 import","import_remote","","bash scripts/import_remote.sh --directory remote/returned/tg2_family --job remote/bundles/raptgen_tg2_family/job.json","Hosted return files -> hash-checked artifacts; no production return exists."),
        ("13 ablation","ablation","ablation","stage ablation scripts/13_structure_ablation.py","Sequence and structure features -> matched held-out scores; test incremental information."),
        ("12 generation / scores","candidates","generation","stage candidates scripts/12_score_candidates.py","Local models -> fixed-budget Markov candidates, novelty and merged local scores; neural joins remain absent."),
        ("14 tertiary","tertiary","","stage tertiary scripts/14_optional_tertiary.py","Optional settings -> explicit exclusion record; no tertiary result is evaluated."),
        ("15 evaluation","evaluation","evaluation","stage evaluation scripts/15_evaluate.py","Fixed test scores -> metrics and family-bootstrap intervals; separate discrimination from sequence coverage."),
        ("16 figures","figures","","stage figures scripts/16_figures.py","Retained tables -> PNG/PDF figures; show all measured arms and unavailable neural branches."),
        ("17 report","report","","stage report scripts/write_report.py","Measured tables and prose template -> README, Markdown and HTML reports; trace numerical claims."),
        ("18 final audit","verification","","stage verification scripts/verify_repository.py --clean-clone-status passed","Verified workflow records -> final audit; only this subprocess is measured here. Run the complete scripts/18_verify.sh first.")]
    stage_records=[]
    for label,stage,state_name,command,purpose in stage_definitions:
        matches=[r for r in resources if r["stage"]==stage]
        basis="latest_recorded_attempt; may_reuse_cache"
        state_path=ROOT/f"logs/state/{state_name}.json"
        if state_name and state_path.exists():
            updated=datetime.fromisoformat(json.loads(state_path.read_text())["time"])
            actual=[r for r in matches if r["exit_status"]=="0" and datetime.fromisoformat(r["start_timestamp"])<=updated<=datetime.fromisoformat(r["end_timestamp"])]
            if actual:matches=actual;basis="successful_execution_containing_scientific_manifest_update"
        r=matches[-1] if matches else None
        stage_records.append(dict(step=label,stage=stage,command=command,input_output_and_reason=purpose,
                                  start_timestamp=r["start_timestamp"] if r else "",elapsed_seconds=r["elapsed_seconds"] if r else "not_run",
                                  peak_rss_bytes=r["peak_rss_bytes"] if r and int(r["peak_rss_bytes"]) else "unavailable",
                                  peak_observed_disk_bytes=r["peak_observed_disk_bytes"] if r else "not_run",exit_status=r["exit_status"] if r else "not_run",
                                  selection=basis if r else "no_production_execution"))
    for r in resources:
        if r["stage"].startswith("preprocess_"):
            stage_records.append(dict(step="05 per-round processing",stage=r["stage"],command="invoked serially by scripts/04_fetch_selex.sh",input_output_and_reason="FASTQ -> verified compact round counts; provider identity and read depth checked.",start_timestamp=r["start_timestamp"],elapsed_seconds=r["elapsed_seconds"],peak_rss_bytes=r["peak_rss_bytes"],peak_observed_disk_bytes=r["peak_observed_disk_bytes"],exit_status=r["exit_status"],selection="retained_per_round_attempt"))
    table(ROOT/"results/pipeline_stage_summary.tsv",stage_records)
    def stage_resource(record):
        if record["elapsed_seconds"]=="not_run":return "Not run"
        mem=f'{int(record["peak_rss_bytes"])/1e6:.1f}' if record["peak_rss_bytes"]!="unavailable" else "unavailable"
        return f'{record["elapsed_seconds"]} s / {mem} MB / {int(record["peak_observed_disk_bytes"])/1e9:.3f} GB'
    stage_table=markdown_table(["Stage","Inputs, outputs and decision","Command","Recorded time / RSS / disk"],[[r["step"],r["input_output_and_reason"],"`"+r["command"]+"`",stage_resource(r)] for r in stage_records[:len(stage_definitions)]])
    novelty=list(rows(ROOT/"results/novelty_distribution.tsv"))
    median_rows={t:[r for r in novelty if r["target"]==t and r["quantity"]=="maximum_training_identity" and float(r["quantile"])==0.5] for t in st}
    def median_range(target):
        numbers=[float(r["value"]) for r in median_rows[target]]
        return f'{min(numbers):.3f}' if min(numbers)==max(numbers) else f'{min(numbers):.3f}-{max(numbers):.3f}'
    def formatted_interval(record):
        return f'{float(record["lower_95"]):+.4f} to {float(record["upper_95"]):+.4f}'
    values=dict(total_candidates=f'{total_candidates:,}',generation_budget=f'{cfg["generation_budget"]:,}',generation_runs=len(generation),
                family_threshold=f'{primary:.0%}',integrin_structure_delta=f'{float(delta("integrin")["difference"]):+.4f}',
                integrin_observed_ap=f'{metric("integrin","enrichment_observed"):.3f}',integrin_ridge_ap=f'{metric("integrin","kmer_ridge"):.3f}',
                tg2_negatives=st["tg2"]["test_negatives"],tg2_test_count=st["tg2"]["test_sequences"],tg2_ridge_ap=f'{metric("tg2","kmer_ridge"):.3f}',tg2_prevalence=f'{metric("tg2","frequency"):.3f}',integrin_prevalence=f'{metric("integrin","frequency"):.3f}',
                sensitivity_thresholds=" and ".join(f'{v:.0%}' for v in cfg["threshold_sensitivity"] if v!=primary),
                data_table=data_table,tg2_nominal_length=design["tg2"]["nominal_length"],integrin_nominal_length=design["integrin"]["nominal_length"],
                tg2_length_interval=design["tg2"]["minimum_length"]+"-"+design["tg2"]["maximum_length"],integrin_length_interval=design["integrin"]["minimum_length"]+"-"+design["integrin"]["maximum_length"],
                tg2_late_unique=f'{int(late_qc["tg2"]["unique_sequences"]):,}',integrin_late_unique=f'{int(late_qc["integrin"]["unique_sequences"]):,}',
                tg2_round_range=design["tg2"]["rounds"].split(",")[0]+"-"+design["tg2"]["rounds"].split(",")[-1],integrin_round_range=design["integrin"]["rounds"].split(",")[0]+"-"+design["integrin"]["rounds"].split(",")[-1],
                assay_total=sum(int(r["retained"]) for r in truth_audit.values()),retrieval_date=retrieval_date,training_cap=f'{cfg["subsampling"]["training_pool"]:,}',
                tg2_pool=f'{int(st["tg2"]["selected_pool"]):,}',integrin_pool=f'{int(st["integrin"]["selected_pool"]):,}',tg2_cap_exclusions=st["tg2"]["pool_cap_exclusions"],integrin_cap_exclusions=st["integrin"]["pool_cap_exclusions"],
                tg2_family_exclusions=st["tg2"]["family_exclusions"],integrin_family_exclusions=st["integrin"]["family_exclusions"],test_table=test_table,ap_table=ap_table,
                bootstrap_replicates=f'{cfg["bootstrap_replicates"]:,}',tg2_bootstrap_valid=f'{int(primary_ci["tg2"]["bootstrap_valid"]):,}',tg2_bootstrap_undefined=primary_ci["tg2"]["bootstrap_undefined"],integrin_bootstrap_valid=f'{int(primary_ci["integrin"]["bootstrap_valid"]):,}',tg2_test_families=st["tg2"]["test_families"],
                split_table=split_table,tg2_split_delta=f'{float(split_delta["tg2"]["difference"]):+.4f}',tg2_split_interval=formatted_interval(split_delta["tg2"]),integrin_split_delta=f'{float(split_delta["integrin"]["difference"]):+.4f}',integrin_split_interval=formatted_interval(split_delta["integrin"]),
                generation_seeds=", ".join(map(str,cfg["generation_seeds"])),curve_budgets=", ".join(f'{n:,}' for n in cfg["generation_curve_budgets"]),
                viennarna_version=system["packages"]["ViennaRNA"],folding_temperature=f'{cfg["folding_temperature_c"]:g}',tg2_structure_delta=f'{float(delta("tg2")["difference"]):+.4f}',tg2_structure_interval=formatted_interval(delta("tg2")),integrin_structure_interval=formatted_interval(delta("integrin")),
                peak_project_gb=f'{peak/1e9:.3f}',disk_budget_gb="13",peak_rss_mb=f'{rss/1e6:.1f}',minimum_free_gb=f'{minfree/1e9:.3f}',resource_table=resource_table,
                python_version=system["packages"]["Python"],bash_version=system["packages"]["Bash"].split("(")[0],
                verification_sentence="Verification is performed with `bash scripts/18_verify.sh`. [Audit records](results/verification.tsv) distinguish local checks from the documented full-mode blocker.")
    values.update(stage_table=stage_table,tg2_identity_medians=median_range("tg2"),integrin_identity_medians=median_range("integrin"),
                  licence_checked_dates=", ".join(sorted({r["licence_checked_date"] for r in rows(ROOT/"results/software_manifest.tsv")} | {r["checked_date"] for r in rows(ROOT/"results/data_terms.tsv")})))
    for t,records in median_rows.items():
        for r in records:
            sources.append(dict(claim_id=f'{t}_seed_{r["seed"]}_median_training_identity',value=r["value"],source_file="results/novelty_distribution.tsv",source_row=f'target={t};method={r["method"]};seed={r["seed"]};quantity={r["quantity"]};quantile={r["quantile"]}',source_column="value",scope="full_candidate_distribution"))
    for r in stage_records:
        if r["elapsed_seconds"]=="not_run":continue
        for col in ["elapsed_seconds","peak_rss_bytes","peak_observed_disk_bytes"]:
            if r[col]=="unavailable":continue
            sources.append(dict(claim_id=f'pipeline_{r["stage"]}_{r["start_timestamp"]}_{col}',value=r[col],source_file="logs/resource_usage.tsv",source_row=f'stage={r["stage"]};start_timestamp={r["start_timestamp"]}',source_column=col,scope=r["selection"]+"; decimal_MB_GB_for_prose"))
    test_summary=ROOT/"results/test_summary.tsv"
    if test_summary.exists():
        test_record=next(rows(test_summary))
        values["verification_sentence"]=f'The recorded test suite has {test_record["passed"]} passing tests, {test_record["failed"]} failures and {test_record["skipped"]} skipped tests. Clean-clone smoke and test checks, shellcheck, restart validation, numerical traceability and leakage checks are recorded in the [verification audit](results/verification.tsv). Full-mode refusal is documented separately from passing local checks.'
        for col in ["passed","failed","skipped"]:
            sources.append(dict(claim_id="tests_"+col,value=test_record[col],source_file="results/test_summary.tsv",source_row="scope="+test_record["scope"],source_column=col,scope="executed_tests"))
    for t,row in design.items():
        for col in ["nominal_length","minimum_length","maximum_length","rounds"]:
            sources.append(dict(claim_id=t+"_design_"+col,value=row[col],source_file="results/dataset_design.tsv",source_row="target="+t,source_column=col,scope="primary_source_configuration"))
        r=late_qc[t]
        sources.append(dict(claim_id=t+"_late_unique_sequences",value=r["unique_sequences"],source_file="results/preprocessing_counts.tsv",source_row=f'target={t};round={r["round"]}',source_column="unique_sequences",scope="measured"))
    for t,r in primary_ci.items():
        for col in ["bootstrap_valid","bootstrap_undefined"]:
            sources.append(dict(claim_id=t+"_primary_"+col,value=r[col],source_file="results/bootstrap_intervals.tsv",source_row=f'target={t};strategy=family;method=kmer_ridge;threshold={primary}',source_column=col,scope="measured"))
    for col in ["retrieval_date_taipei","python_version","bash_version","viennarna_version"]:
        sources.append(dict(claim_id=col,value=provenance[col],source_file="results/readme_provenance.tsv",source_row="scope=readme_sources",source_column=col,scope="manifest_or_installed_version"))
    for r in executed:
        for col in ["elapsed_seconds","peak_rss_bytes","peak_observed_disk_bytes"]:
            sources.append(dict(claim_id="scientific_execution_"+r["stage"]+"_"+col,value=r[col],source_file="logs/resource_usage.tsv",source_row=f'stage={r["stage"]};start_timestamp={r["start_timestamp"]}',source_column=col,scope="successful_uncached_execution; decimal_MB_GB_for_prose"))
    table(ROOT/"results/readme_traceability.tsv",sources)
    readme=render_template((ROOT/"docs/readme_template.md").read_text(encoding="utf-8"),values)

    atomic_text(ROOT/"README.md",readme)
    report="# Aptamer SELEX analysis report\n"+readme.partition("\n")[2]
    report=re.sub(r"(!?\[[^\]]*\]\()([^)]+)(\))",lambda m:m[1]+(m[2] if re.match(r"https?://|#",m[2]) else "../"+m[2])+m[3],report)
    atomic_text(ROOT/"docs/analysis_report.md",report)
    atomic_text(ROOT/"results/report.html",report_html(readme))
    atomic_text(ROOT/"HANDOVER.md",'''# Remaining external work

The public sequencing processing, primary assay reconstruction, local baselines, leakage audits, Markov generation, secondary-structure ablation, statistics and figures are real outputs. No RaptGen or AptaDiff model has trained in this run. No neural ranking, neural recovery or tertiary result exists.

RaptGen needs a hosted Linux machine with a GPU compatible with the pinned legacy stack. Its environment resolution and training must be verified there. Prepare with `bash run_all.sh --mode core --from 10`. On the host run `conda env create -f remote/raptgen/environment.yml`, activate it, then `bash remote/raptgen/run_remote.sh remote/bundles/raptgen_tg2_family output/tg2_family`. Repeat both targets and the matched diagnostic. Return only the files listed in remote/README.md and import with `bash scripts/import_remote.sh --directory remote/returned/tg2_family --job remote/bundles/raptgen_tg2_family/job.json`. Source/configuration/input/output hashes are independently checked. The deleted checkpoint's hash is remote-recorded rather than independently recomputed locally.

AptaDiff needs explicit software and weight licensing, then an authorized implementation of its published VAE/conditional-diffusion interface. Its present environment is a preflight placeholder and its launcher refuses execution. See remote/aptadiff/BLOCKER.md. No substitute model is presented as AptaDiff. InstructNA remains excluded for unclear licensing and unverified weight/resource terms.

The first exact-only diagnostic was mislabeled as an independent random split. The final repository retains it as a matched exclusion diagnostic and adds a separate per-sequence random partition. This correction is declared in configuration and methods. Primary test IDs and generation seeds were preserved. Initial installation and metadata RSS rows were unavailable under the first wrapper; their original logs remain. Later production stages use the corrected sampler. Direct DDBJ metadata routes failed; DDBJ-submitted INSDC metadata was verified through ENA. Those failures are recorded.

The local artifact importer is implemented and tested, but stage 12 currently merges local scores only. After the first real neural artifacts return, implement and validate their joins to local assay records, candidate recovery/novelty aggregation, seed statistics and report coverage. The full preflight explicitly refuses this missing integration even if returned files appear, so it cannot silently complete with only baselines. This part has not been executed or validated against actual hosted outputs.

After the neural licensing/implementation, returned artifacts and aggregation are resolved, resume with `bash run_all.sh --mode full --gpu-mode hosted --from 11`, then run `bash scripts/18_verify.sh`. Full mode currently stops explicitly. The repository is an executed local core benchmark with a documented full-benchmark blocker, not a completed neural comparison or a finished publication claim.
''')
    metadata_path=ROOT/"docs/github_metadata.json"
    metadata=json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    metadata.update(description="Family-held-out aptamer benchmark on two public HT-SELEX datasets, with experimental discrimination, a Markov generation baseline, and canonical-RNA structure ablation.",topics=["aptamer","ht-selex","rna","sequence-analysis","secondary-structure","viennarna","machine-learning","benchmark","bioinformatics","computational-biology","nucleic-acids","selex","reproducibility","python","bash"])
    metadata.setdefault("status","prepared_local_metadata_no_remote_configured")
    atomic_text(metadata_path,json.dumps(metadata,indent=2))
    print("README and analysis report written from measured tables")


if __name__=="__main__":
    main()
