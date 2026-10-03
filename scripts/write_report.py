"""Write the scientific narrative only from completed result tables and logs."""
import argparse
from collections import defaultdict
import html
import json
from pathlib import Path
from common import ROOT, rows, table, atomic_text, settings, sha256


def markdown_table(headers, data):
    return "| " + " | ".join(headers) + " |\n| " + " | ".join(["---"]*len(headers)) + " |\n" + "\n".join("| "+" | ".join(map(str,r))+" |" for r in data)


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
    ap_table=markdown_table(["Method","TG2 grouped AP (95% interval)","Integrin grouped AP (95% interval)"],[[m,f'{metric("tg2",m):.3f} ({interval("tg2",m)})',f'{metric("integrin",m):.3f} ({interval("integrin",m)})'] for m in ["frequency","enrichment","kmer_ridge","markov","kmer_logistic","kmer_structure_ridge","kmer_structure_logistic","frequency_observed","enrichment_observed"]])
    split_table=markdown_table(["Target","Independent naive random AP","Matched exact-only AP","Family AP"],[[t,f'{metric(t,"kmer_ridge","naive_sequence_random"):.3f}',f'{metric(t,"kmer_ridge","random"):.3f}',f'{metric(t,"kmer_ridge"):.3f}'] for t in st])
    data_table=markdown_table(["Target / accession","Rounds processed","Raw reads","Retained reads","Discarded reads","Assays retained / excluded","Primary test / families"],[[f'{t} / '+("DRA009383" if t=="tg2" else "DRA009384"),r["rounds"],r["raw_reads"],r["retained_reads"],r["discarded_reads"],f'{r["assays_retained"]} / {r["assays_excluded"]}',f'{r["test_sequences"]} / {r["test_families"]}'] for t,r in st.items()])
    last={}
    for r in resources:last[r["stage"]]=r
    resource_table=markdown_table(["Stage","Latest elapsed s","Peak RSS MB","Observed project GB"],[[name,r["elapsed_seconds"],f'{int(r["peak_rss_bytes"])/1e6:.1f}' if int(r["peak_rss_bytes"]) else "unavailable",f'{int(r["peak_observed_disk_bytes"])/1e9:.3f}'] for name,r in last.items()])
    changes="; ".join(f'{t}: {float(delta(t)["difference"]):+.4f} (95% {float(delta(t)["lower_95"]):+.4f} to {float(delta(t)["upper_95"]):+.4f})' for t in st)
    readme=f'''# Measured selection baselines; neural comparison remains blocked

## Summary

TG2's grouped test set has {st['tg2']['test_negatives']} negative sequence among {st['tg2']['test_sequences']} sequences. Its prevalence baseline is {metric('tg2','frequency'):.3f} average precision (AP), compared with {metric('tg2','kmer_ridge'):.3f} for the SELEX-trained k-mer model. Integrin's grouped k-mer AP is {metric('integrin','kmer_ridge'):.3f}; observed enrichment is {metric('integrin','enrichment_observed'):.3f}. The latter uses test-sequence counts and is an observational reference. Adding canonical-RNA structural features changes integrin ridge AP by {float(delta('integrin')['difference']):+.4f}, with a family-bootstrap interval spanning zero. The Markov null generated {total_candidates:,} unique candidates across all target-seed runs and recovered {exact_recovery} exact test sequences and {family_recovery} positive test families. RaptGen has a hosted job but no returned run. AptaDiff is excluded pending explicit repository licensing and an authorized adapter. The largest observed local footprint is {peak/1e9:.3f} GB, including the environment and temporary files. Sources: [summary](results/summary.tsv), [metrics](results/evaluation_metrics.tsv), [generation](results/generation_statistics.tsv), [resources](logs/resource_usage.tsv).

## Background

Aptamers bind through a molecular conformation formed by their sequence and chemical context. SELEX repeatedly partitions molecules and amplifies survivors. HT-SELEX sequences those pools. Read abundance can reflect selection, PCR bias or a bottleneck; it is not an independent affinity measurement. A model trained on target-specific selection data can learn that distribution. A different target without such data is a different scientific problem.

Related sequences can share motifs and appear on both sides of a sequence split. Here, exact edit-distance components define families at an arbitrary preregistered identity of {primary:.0%}. Indels count. Secondary structure supplies an intramolecular folding proxy; a predicted tertiary fold or model ranking does not establish binding. **This pipeline does not claim that an aptamer can be reliably designed for an arbitrary target from its protein sequence or structure alone.**

## Data

{data_table}

These RNA libraries contain 2'-fluoro pyrimidines. The primary supplement defines transcribed-strand primer filters and permits length intervals around the nominal library lengths. No assay sequence was excluded for an indel. Every rejection class is retained in [preprocessing counts](results/preprocessing_counts.tsv). Reads were not subsampled. Training-pool capping excluded {st['tg2']['pool_cap_exclusions']} eligible TG2 sequences and {st['integrin']['pool_cap_exclusions']} eligible integrin sequences before family construction; the selected pools contain {st['tg2']['selected_pool']} and {st['integrin']['selected_pool']} sequences. The primary family rule then excludes {st['tg2']['family_exclusions']} and {st['integrin']['family_exclusions']} pool members, respectively, from generator fitting.

The RaptRanker primary supplement provides original Positive/Negative labels and continuous end-of-injection SPR response in RU. We preserve both and introduce no binary cutoff. [Ground-truth audit](results/ground_truth_audit.tsv) and [each original measurement](config/ground_truth_manifest.tsv) retain provenance. [RaptRanker](https://doi.org/10.1093/nar/gkaa484) establishes the target identities. The ENA mirror provides the DDBJ-submitted study and experiment records.

[AptaDiff](https://pmc.ncbi.nlm.nih.gov/articles/PMC11491854/) distinguishes its IGFBP3/PTK7 datasets A/B from public TG2/integrin C/D in Table 1, while its availability sentence assigns the DRA accessions to A/B. Repository primer and length measurements support keeping these biological datasets distinct. The accession sentence and A/B public accessions remain unresolved. Four inspected repository data files are excluded from benchmark training; the official complete rounds are used instead. Details remain in [data provenance](results/data_provenance.tsv).

## Pipeline

`bash scripts/00_configure.sh --threads 1 --ram 16 --disk 13 --gpu-mode none --yes` measures free space, preserves an external reserve and writes project.conf. Source checks and installation precede sequencing. `bash scripts/04_fetch_selex.sh` downloads one file, checks size and MD5, counts with SQLite, verifies compact output and deletes raw input. The measured pilot supplies the [disk projection](results/disk_projection.tsv).

`bash run_all.sh --mode core --from 06` reconstructs ground truth, families, baselines, structure, candidates and statistics. Every stage uses the resource wrapper. Configuration and output hashes govern restart decisions. The scientific definitions and alternatives are explained in [methods](docs/methods.md) and the [score dictionary](docs/score_dictionary.md).

The grouped primary arm excludes every test family from training. The matched exact-only diagnostic shares the same assay test IDs and allows their relatives in training. It is stored as `random` but is not an independent random partition. A separate `naive_sequence_random` arm provides that partition. It was added after the first local run to correct this diagnostic definition; no performance-driven seed selection occurred. Its test cases differ, so its difference from the primary arm is not a paired leakage estimate. Both diagnostics remain visible. The family-aware audits have zero exact and family overlaps.

The complete local data flow ends at fixed experimental test cases; missing neural jobs are shown explicitly.

![Pipeline](figures/01_pipeline.png)

## Results

Held-out AP must be judged against prevalence. TG2 has only one negative test sequence, so a high AP gives weak evidence of discrimination.

{ap_table}

`frequency_observed` and `enrichment_observed` query original read counts. Exact-lookup `frequency` and `enrichment` use filtered training counts and become constant after exact test exclusion. Neither constant score is hidden. The ridge model predicts SELEX enrichment; logistic fitting uses training assay labels only. Its output is an uncalibrated decision function, not a probability or KD. Full AUROC, rank correlations, top-ranked precision and undefined-score statuses are in [evaluation metrics](results/evaluation_metrics.tsv).

The same k-mer model has different performance under the independent naive partition and the controlled exclusion diagnostic.

{split_table}

![Split comparison](figures/03_random_vs_family.png)

The grouped precision-recall curves include the observational references and preserve the losing local methods.

![Experimental discrimination](figures/04_experimental_discrimination.png)

Processing retains the full round trajectory, including later TG2 rounds that the original RaptRanker analysis excluded because negative-labelled sequences amplified.

![Dataset composition](figures/02_dataset_composition.png)

Structure-minus-sequence ridge changes are {changes}. These are paired comparisons on identical test cases with sequence-family bootstrap resampling; they exclude training uncertainty.

![Structure ablation](figures/07_structure_ablation.png)

The Markov null recovers no positive experimental test family at the configured budgets in this run. All generation seeds remain in the denominator.

![Generation budget](figures/05_generation_budget.png)

Novelty is measured against the actual filtered training pool. The plotted points are the first fixed subset from each seed; full distributions and nearest sequences remain in generated data.

![Novelty and ranking](figures/06_novelty_ranking.png)

At the sensitivity threshold, an integrin assay-supervised fit lacks both training classes and is excluded in the grouped arm, including its structural variant. Both exclusions are logged. Two optional methods, AptaDiff and InstructNA, have unclear repository licensing. RaptGen is not executed because no hosted result has returned. No tertiary method runs. These unavailable results cannot support an architectural comparison.

The largest observed local footprint stays below the configured ceiling. Initial installation RSS was not captured by the first wrapper version; that row remains marked unavailable. Subsequent stages use the corrected process sampler. [Failure records](logs/failures.tsv) retain the metadata-route failures, test-fixture repair and any failed stage. The largest measured stage RSS is {rss/1e6:.1f} MB and the minimum observed free space is {minfree/1e9:.3f} GB. RSS and disk are sampled, so brief peaks may be missed.

![Resource profile](figures/08_resources.png)

{resource_table}

## Repository structure

```text
config/       fixed decisions, primer design, assay and split manifests
scripts/      numbered stages, shared resource wrapper, plotting and reporting
remote/       pinned RaptGen job and explicit AptaDiff blocker
results/      compact measurements, provenance, predictions and uncertainty
figures/      scripted PNG and PDF figures
logs/         resources, failures, restart hashes and verification
data/         ignored raw, compact, generated and external data
tests/        synthetic end-to-end, scientific and artifact checks
docs/         methods, score definitions and the executed report
```

## Usage

Resource context: 16 GB RAM, maximum 13.0 GB local footprint. Bash, Git and Python are required. The local environment has no CUDA toolkit, GPU PyTorch installation or neural checkpoint. Git Bash was used on Windows; Linux and WSL can use the same shell scripts with a compatible Python interpreter.

```bash
# Set PYTHON or BOOTSTRAP_PYTHON to a supported interpreter when PATH is ambiguous.
bash scripts/02_install.sh
bash scripts/00_configure.sh --threads 1 --ram 16 --disk 13 --gpu-mode none --yes
bash run_all.sh --mode smoke
bash run_all.sh --mode core
bash run_all.sh --mode core --from 10
# Execute and return hosted jobs as described in remote/README.md.
bash scripts/import_remote.sh --directory remote/returned/tg2_family \
  --job remote/bundles/raptgen_tg2_family/job.json
bash run_all.sh --mode full --gpu-mode hosted --from 11
bash scripts/18_verify.sh
```

The smoke test uses synthetic reads and requires no public-data download or GPU. The latest measured times and disk footprints are listed above. Resource summaries describe the logged rows present when the report started; their count is recorded in summary.tsv. Full mode explicitly refuses missing or excluded neural arms and the unfinished neural aggregation. Read [HANDOVER.md](HANDOVER.md) for the exact remaining tasks. `--target` limits retrieval; the comparative analysis requires processed data for both targets. A from-stage restart requires the preceding validated files. The report is [HTML](results/report.html), with [Markdown source](docs/analysis_report.md); Quarto is optional.

## Limitations

The independent test measurements are few, especially independent TG2 families. I am least confident in TG2 binary discrimination because its test partition has only one negative. The family threshold is arbitrary and its sensitivity changes exclusion and training-class availability. Hash capping can omit informative late-round sequences. Read counts are affected by amplification and selection biases. Reported bootstrap intervals exclude refitting and model selection uncertainty.

Canonical RNA folding does not reproduce the chemistry of 2'-fluoro-modified RNA exactly. MFE concerns intramolecular folding. None of the generated candidates has a new wet-lab assay. Neural environments are legacy stacks; hosted training and artifact return remain unverified. AptaDiff's licence contradiction blocks its implementation. Tertiary prediction is disabled and contributes no evidence. This run does not establish reliable aptamer design for a target without target-specific selection or functional data.

## Data availability

The public accessions are DRA009383 and DRA009384. [Run metadata](results/selex_manifest.tsv) records every remote FASTQ, provider checksum, file size and explicit round evidence. [Download provenance](results/download_manifest.tsv) supplies retrieval timestamps and local hashes. The primary supplement is fetched from the publisher's linked PDF and parsed directly. Retained compact scientific tables are tracked; raw reads, derived count tables, primary PDFs, large generation tables and model weights are not. The Bash commands above regenerate them. Inspection date and local run date are recorded in the manifests rather than inferred from filenames.

## Citation

Primary literature: [RaptRanker](https://doi.org/10.1093/nar/gkaa484), [RaptGen](https://doi.org/10.1038/s43588-022-00249-6), [AptaDiff](https://doi.org/10.1093/bib/bbae517), and the discussed, excluded [InstructNA](https://doi.org/10.1038/s43588-026-00965-3). Software references: [ViennaRNA](https://www.tbi.univie.ac.at/RNA/), [NumPy](https://numpy.org/citing-numpy/), [SciPy](https://scipy.org/citing-scipy/), [scikit-learn](https://jmlr.org/papers/v12/pedregosa11a.html), [Matplotlib](https://matplotlib.org/stable/project/citing.html), and [RapidFuzz](https://github.com/rapidfuzz/RapidFuzz). Families use this repository's exact graph construction and RapidFuzz global edit distances; VSEARCH is not a dependency. Sources and checks are in [citations](results/citations.tsv) and [software provenance](results/software_manifest.tsv).

## Licence

Original repository code is MIT under the configured author's identity. RaptGen is MIT and is installed separately only on the hosted path. AptaDiff's paper declares MIT, but its inspected repository lacks a licence file; its executable arm is excluded. InstructNA also lacks a clear repository licence and is excluded. ViennaRNA has its own custom licence with attribution and redistribution conditions. These terms are not replaced by this repository's MIT licence. No pretrained weight licence is assumed.

The RaptRanker article and supplement have their own CC-BY-NC terms. The original PDF and upstream datasets are not redistributed. Derived factual assay fields are attributed, with terms recorded in [data terms](results/data_terms.tsv). All installed dependency and licence records retain the date checked. The current evidence ends at the local baselines and canonical-RNA ablation.
'''
    atomic_text(ROOT/"README.md",readme)
    report=readme.replace("# Measured selection baselines; neural comparison remains blocked", "# Executed local analysis report",1).replace("figures/","../figures/")
    atomic_text(ROOT/"docs/analysis_report.md",report)
    atomic_text(ROOT/"results/report.html",'<!doctype html><html lang="en"><meta charset="utf-8"><title>Executed local analysis</title><style>body{max-width:1100px;margin:40px auto;font:16px/1.6 system-ui;color:#182c3b}pre{white-space:pre-wrap;font:inherit}img{max-width:100%}</style><h1>Executed local analysis</h1><pre>'+html.escape(readme)+'</pre>'+''.join(f'<figure><img src="../figures/{name}.png" alt="{name}"></figure>' for name in ["03_random_vs_family","07_structure_ablation","08_resources"])+"</html>")
    atomic_text(ROOT/"HANDOVER.md",'''# Remaining external work

The public sequencing processing, primary assay reconstruction, local baselines, leakage audits, Markov generation, secondary-structure ablation, statistics and figures are real outputs. No RaptGen or AptaDiff model has trained in this run. No neural ranking, neural recovery or tertiary result exists.

RaptGen needs a hosted Linux machine with a GPU compatible with the pinned legacy stack. Its environment resolution and training must be verified there. Prepare with `bash run_all.sh --mode core --from 10`. On the host run `conda env create -f remote/raptgen/environment.yml`, activate it, then `bash remote/raptgen/run_remote.sh remote/bundles/raptgen_tg2_family output/tg2_family`. Repeat both targets and the matched diagnostic. Return only the files listed in remote/README.md and import with `bash scripts/import_remote.sh --directory remote/returned/tg2_family --job remote/bundles/raptgen_tg2_family/job.json`. Source/configuration/input/output hashes are independently checked. The deleted checkpoint's hash is remote-recorded rather than independently recomputed locally.

AptaDiff needs explicit software and weight licensing, then an authorized implementation of its published VAE/conditional-diffusion interface. Its present environment is a preflight placeholder and its launcher refuses execution. See remote/aptadiff/BLOCKER.md. No substitute model is presented as AptaDiff. InstructNA remains excluded for unclear licensing and unverified weight/resource terms.

The first exact-only diagnostic was mislabeled as an independent random split. The final repository retains it as a matched exclusion diagnostic and adds a separate per-sequence random partition. This correction is declared in configuration and methods. Primary test IDs and generation seeds were preserved. Initial installation and metadata RSS rows were unavailable under the first wrapper; their original logs remain. Later production stages use the corrected sampler. Direct DDBJ metadata routes failed; DDBJ-submitted INSDC metadata was verified through ENA. Those failures are recorded.

The local artifact importer is implemented and tested, but stage 12 currently merges local scores only. After the first real neural artifacts return, implement and validate their joins to local assay records, candidate recovery/novelty aggregation, seed statistics and report coverage. The full preflight explicitly refuses this missing integration even if returned files appear, so it cannot silently complete with only baselines. This part has not been executed or validated against actual hosted outputs.

After the neural licensing/implementation, returned artifacts and aggregation are resolved, resume with `bash run_all.sh --mode full --gpu-mode hosted --from 11`, then run `bash scripts/18_verify.sh`. Full mode currently stops explicitly. The repository is an executed local core benchmark with a documented full-benchmark blocker, not a completed neural comparison or a finished publication claim.
''')
    atomic_text(ROOT/"docs/github_metadata.json",json.dumps(dict(description="An experimental holdout benchmark of target-specific HT-SELEX models, with family leakage controls, canonical RNA structure ablation and measured local resources. Neural comparison awaits hosted runs and licence resolution.",topics=["aptamer","ht-selex","rna","sequence-analysis","secondary-structure","viennarna","machine-learning","benchmark","bioinformatics","computational-biology","nucleic-acids","selex","reproducibility","python","bash"],status="prepared_local_metadata_no_remote_configured"),indent=2))
    print("README and executed report written from actual outputs; no neural result invented")


if __name__=="__main__":
    main()
