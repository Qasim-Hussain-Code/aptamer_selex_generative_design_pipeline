# No positive test-family recovery at 10,000 Markov candidates per run

## Summary

This benchmark tests sequence ranking and generation against published SPR assay labels for TG2 and integrin alpha V beta 3, holding out families at 80% edit identity. SELEX enrichment-ridge average precision (AP) is 0.926 and 0.744, against constant-frequency baselines of 0.933 and 0.600. Integrin observed enrichment reaches 0.875 using test-sequence read counts. Matched exact-only exclusion changes ridge AP by -0.0143 for TG2 and +0.0457 for integrin; both paired intervals include zero. Across 10 target-seed runs, the Markov model produces 100,000 candidates with 10,000 unique sequences per run, recovering no exact test sequence or positive test family. Structure changes integrin ridge AP by +0.0031, with an interval spanning zero. Peak sampled disk use is 0.546 GB. TG2 has only 1 negative among 15 test sequences. RaptGen training is unrun; AptaDiff licensing and neural aggregation remain unresolved. [Summary tables](results/summary.tsv), [evaluation metrics](results/evaluation_metrics.tsv) and [paired comparisons](results/paired_differences.tsv) contain the measurements.

## Background

Aptamers are nucleic-acid molecules whose sequence-dependent folds can recognize a target. SELEX repeatedly selects a library against that target and amplifies the retained molecules; HT-SELEX adds sequencing across rounds. Read abundance reflects selection, amplification bias and sampling, so enrichment can differ from experimentally measured affinity. A variable region of length L has 4^L possible nucleotide sequences, far more than the 6,031 and 12,250 distinct final-round sequences observed here.

A selection-trained generator learns the sampled sequence distribution and its biases for a particular target and protocol. Its probability score does not validate binding. This benchmark asks whether that information ranks independently assayed sequences after related training sequences have been excluded. Close relatives across a split can make sequence memorization appear to generalize; holding out families tests a more demanding form of transfer within these target-specific datasets.

Secondary structure constrains which nucleotides are paired and which remain available for recognition. The ablation tests whether those descriptors add information to sequence features. Tertiary prediction requires long-range packing and target-interaction assumptions that are less established for modified aptamers; a predicted fold alone supplies no functional validation. No tertiary model runs here, and the canonical-RNA approximation does not reproduce the libraries' modified chemistry exactly.

Families are connected components of global edit identity, with insertions and deletions included. The primary threshold of 80% was fixed before model evaluation; 70% and 90% are retained as sensitivity analyses. The threshold is arbitrary. This pipeline does not claim that an aptamer can be reliably designed for an arbitrary target from its protein sequence or structure alone.

## Data

| Target / accession | Rounds | Raw reads | Retained reads | Rejected reads | Assays retained / excluded |
| --- | --- | --- | --- | --- | --- |
| TG2 / DRA009383 | 0-8 | 840,620 | 650,294 | 190,326 | 30 / 0 |
| Integrin alpha V beta 3 / DRA009384 | 3-6 | 467,095 | 306,561 | 160,534 | 48 / 0 |

TG2 variable regions have a nominal length of 30 nt and a retained interval of 25-35 nt; integrin regions have a nominal length of 40 nt and a retained interval of 35-45 nt. Both libraries contain 2'-fluoro pyrimidines. All public reads were processed, with no read subsampling. The final processed rounds contain 6,031 unique TG2 sequences and 12,250 unique integrin sequences. Counts and rejection classes for each round are recorded in [preprocessing QC](results/preprocessing_counts.tsv).

After primer and length filtering, the final rounds retain 6,031 distinct TG2 sequences and 12,250 integrin sequences ([QC counts](results/preprocessing_counts.tsv)).

![Read counts and unique sequences across the TG2 and integrin SELEX rounds](figures/02_dataset_composition.png)

Figure 2. Raw reads, retained reads and unique retained sequences for every processed round. TG2 rounds 0-8 and integrin rounds 3-6 are retained, including late TG2 rounds omitted from the original RaptRanker ranking analysis.

The 78 assay records come directly from RaptRanker Supplementary Tables S4/S5. Original Positive/Negative labels and end-of-injection SPR responses in response units (RU) are preserved. No new binary cutoff was introduced, and no assay was excluded for an indel. The [ground-truth manifest](config/ground_truth_manifest.tsv) records each source measurement and its provenance. Public reads and reference inputs were retrieved on 2026-10-03 (Asia/Taipei); exact URLs, checksums and retrieval timestamps are in the [download manifest](results/download_manifest.tsv).

AptaDiff Table 1 identifies its datasets A/B as IGFBP3/PTK7 and its public datasets C/D as TG2/integrin, while its data-availability sentence associates the DRA accessions with A/B. The accession records and repository primer checks support the TG2/integrin mapping used here. Ambiguous repository inputs are excluded; the original contradiction remains in [data provenance](results/data_provenance.tsv).

There is also a literal integrin primer-boundary discrepancy: the article's PCR forward primer ends in `CAGAA`, while Supplementary Table S3's computational filter ends in `CAGAAG`. Processing follows S3, and assay strings follow S5. Exact count lookup and sequence recovery depend on this convention. Primer definitions and the rationale are recorded in [dataset configuration](config/datasets.tsv) and [methods](docs/methods.md).

## Pipeline

The workflow installs the pinned CPU environment, verifies primary sources, resolves sequencing rounds from experiment metadata, then downloads and processes one FASTQ file at a time. Each file is checked against the provider's size, MD5 and read count. Retained variable regions are counted with SQLite, compact outputs are reconciled against QC totals, and the raw file is removed. Hashes of configuration, scientific code, inputs and outputs control analysis restarts.

The executed workflow evaluates local SELEX fits and canonical-RNA ablations, while neural branches remain without returned results ([figure manifest](results/figure_manifest.tsv)).

![Workflow from public sequencing reads through family exclusion, baselines, structure and experimental evaluation](figures/01_pipeline.png)

Figure 1. Data processing and evaluation workflow. Experimental responses are isolated from SELEX-only fits and generation. The neural branches have no measured results in this repository.

Late-round sequences observed at least twice are selected by seeded hash, with a cap of 2,500 sequences per target. This retains 1,893 TG2 sequences and 2,500 integrin sequences; the cap excludes 0 and 59 eligible sequences, respectively. The family-aware primary split then removes 83 TG2 and 298 integrin pool members from fitting. Within the selected pool and assay-sequence graph, both primary targets have zero exact and family overlaps between training and test.

Three split definitions are recorded in the [split manifest](config/split_manifest.tsv):

- `family` holds out assay families and excludes their relatives from fitting.
- `random` is a matched exact-only exclusion diagnostic using the same test sequences as `family`; it is not an independent random partition.
- `naive_sequence_random` partitions assay sequences independently by seeded hash and permits relatives in training. This arm was added after the first local run to correct the diagnostic definition, without selecting a seed or metric for performance.

The ridge model predicts final-versus-previous-round log2 enrichment. The logistic model instead uses labels from training assays. Both fit feature scaling on their training data only. Hyperparameters, seeds, score definitions and structural descriptors are recorded in [benchmark configuration](config/benchmark.yml), [methods](docs/methods.md) and the [score dictionary](docs/score_dictionary.md).

After installation, `source scripts/lib.sh` selects the project interpreter and defines the resource-measured `stage` command. Individual commands below require the inputs from preceding stages; `run_all.sh` supplies that order and acquires the reference inputs.

| Stage | Inputs, outputs and decision | Command | Recorded time / RSS / disk |
| --- | --- | --- | --- |
| 00 configure | Host resources -> project.conf; enforce memory, disk and thread limits. | `bash scripts/00_configure.sh --threads 1 --ram 16 --disk 13 --gpu-mode none --yes` | 1.408 s / 28.6 MB / 0.485 GB |
| 01 sources | Acquired publications and repository metadata -> provenance and licence records; establish source identity. | `stage sources scripts/01_verify_sources.py --offline` | 1.731 s / 35.5 MB / 0.483 GB |
| 02 install | Pinned package versions -> CPU environment; keep dependency and chemistry assumptions explicit. | `bash scripts/02_install.sh` | 190.033 s / 182.9 MB / 0.546 GB |
| 03 metadata | Accessions and experiment XML -> run manifest; resolve rounds from source evidence. | `stage metadata scripts/03_fetch_metadata.py` | 0.973 s / 6.2 MB / 0.485 GB |
| 04/05 reads | One FASTQ at a time -> verified compact counts and QC; displayed time is cached validation, with per-round processing listed separately. | `bash scripts/04_fetch_selex.sh --target all` | 0.966 s / 4.1 MB / 0.485 GB |
| 06 assays | Supplement S4/S5 -> assay manifest; preserve original labels and SPR units. | `stage ground_truth scripts/06_build_ground_truth.py` | 2.395 s / 53.1 MB / 0.481 GB |
| 07 splits | Late-round pool and assays -> family partitions; exclude test relatives before fitting. | `stage splits scripts/07_build_splits.py` | 9.510 s / 79.6 MB / 0.482 GB |
| 08 baselines | Training counts and permitted assays -> fixed sequence fits and scores; retain exposure-matched baselines. | `stage baselines scripts/08_baselines.py` | 24.918 s / 205.0 MB / 0.481 GB |
| 09 structure | Complete assay constructs -> canonical-RNA descriptors; approximate intramolecular structure. | `stage secondary scripts/09_secondary_structure.py` | 29.962 s / 44.0 MB / 0.481 GB |
| 10 GPU jobs | Training-only exports -> pinned RaptGen job bundles; isolate held-out assay responses. | `stage gpu_jobs scripts/10_prepare_gpu_jobs.py` | 1.383 s / 63.3 MB / 0.485 GB |
| 11 import | Hosted return files -> hash-checked artifacts; no production return exists. | `bash scripts/import_remote.sh --directory remote/returned/tg2_family --job remote/bundles/raptgen_tg2_family/job.json` | Not run |
| 13 ablation | Sequence and structure features -> matched held-out scores; test incremental information. | `stage ablation scripts/13_structure_ablation.py` | 25.520 s / 211.4 MB / 0.481 GB |
| 12 generation / scores | Local models -> fixed-budget Markov candidates, novelty and merged local scores; neural joins remain absent. | `stage candidates scripts/12_score_candidates.py` | 48.516 s / 145.9 MB / 0.481 GB |
| 14 tertiary | Optional settings -> explicit exclusion record; no tertiary result is evaluated. | `stage tertiary scripts/14_optional_tertiary.py` | 0.965 s / 6.2 MB / 0.485 GB |
| 15 evaluation | Fixed test scores -> metrics and family-bootstrap intervals; separate discrimination from sequence coverage. | `stage evaluation scripts/15_evaluate.py` | 25.406 s / 131.3 MB / 0.481 GB |
| 16 figures | Retained tables -> PNG/PDF figures; show all measured arms and unavailable neural branches. | `stage figures scripts/16_figures.py` | 5.256 s / 238.8 MB / 0.485 GB |
| 17 report | Measured result tables -> README, Markdown and HTML reports; trace numerical claims. | `stage report scripts/write_report.py` | 1.160 s / 36.6 MB / 0.484 GB |
| 18 final audit | Verified workflow records -> final audit; only this subprocess is measured here. Run the complete scripts/18_verify.sh first. | `stage verification scripts/verify_repository.py --clean-clone-status passed` | 13.209 s / 85.0 MB / 0.484 GB |

Commands run from the repository root. The [stage execution table](results/pipeline_stage_summary.tsv) joins these inputs, outputs and decisions to measured elapsed time, sampled RSS and sampled project size. It identifies successful scientific executions where a manifest is available; other rows use the latest recorded attempt and label that basis. Downloads are measured per FASTQ, and hosted import has no production execution. Resource rows describe the recorded execution, including cache reuse where declared.

## Results

### Experimental discrimination

| Primary panel | Test sequences | Positive / negative | Independent test families |
| --- | --- | --- | --- |
| TG2 | 15 | 14 / 1 | 6 |
| Integrin alpha V beta 3 | 25 | 15 / 10 | 24 |

| Method | TG2 AP (95% interval) | Integrin AP (95% interval) |
| --- | --- | --- |
| Training frequency lookup | 0.933 (0.692 to 0.957) | 0.600 (0.400 to 0.792) |
| Training enrichment lookup | 0.933 (0.692 to 0.957) | 0.600 (0.400 to 0.792) |
| SELEX enrichment ridge | 0.926 (0.639 to 0.965) | 0.744 (0.526 to 0.971) |
| Markov log probability per nt | 0.870 (0.782 to 0.910) | 0.745 (0.524 to 0.920) |
| Assay-label logistic | 0.957 (0.748 to 0.992) | 0.568 (0.380 to 0.862) |
| SELEX ridge + structure | 0.926 (0.642 to 0.965) | 0.747 (0.533 to 0.974) |
| Assay logistic + structure | 0.957 (0.748 to 0.992) | 0.574 (0.382 to 0.860) |
| Observed frequency | 0.955 (0.825 to 0.982) | 0.765 (0.549 to 0.945) |
| Observed enrichment | 0.957 (0.844 to 0.983) | 0.875 (0.702 to 0.977) |

The exact-lookup frequency and enrichment baselines become constant after test-sequence exclusion, so their AP equals test prevalence. `frequency_observed` and `enrichment_observed` query original read counts, including counts of test sequences. Their exposure differs from that of the held-out models. Ridge scores predict selection enrichment; logistic scores are uncalibrated decision functions. Neither is a KD estimate.

Intervals resample sequence families for 2,000 replicates with one fixed fitted model. TG2 retains 1,314 two-class resamples and excludes 686 single-class resamples; integrin retains 2,000. The intervals are conditional on retaining both classes and exclude refitting uncertainty. TG2 has only 6 independent test families. [Bootstrap records](results/bootstrap_intervals.tsv) preserve the undefined counts; secondary endpoints, including AUROC, response rank correlations and top-ranked precision, remain in [evaluation metrics](results/evaluation_metrics.tsv).

On the primary panels, ridge AP is 0.926 for TG2 and 0.744 for integrin; TG2's 0.933 prevalence exceeds its ridge AP ([metrics](results/evaluation_metrics.tsv)).

![Precision-recall curves for the family-held-out experimental TG2 and integrin panels](figures/04_experimental_discrimination.png)

Figure 4. Precision-recall curves on the primary experimental panels. Observed frequency and enrichment are exposed references. The high positive prevalence of the TG2 panel limits what its AP can establish about discrimination.

### Split comparisons

| Target | Split | Test sequences | Positive fraction | Ridge AP |
| --- | --- | --- | --- | --- |
| TG2 | Independent sequence partition | 15 | 0.800 | 0.739 |
| TG2 | Matched exact-only exclusion | 15 | 0.933 | 0.912 |
| TG2 | Family exclusion | 15 | 0.933 | 0.926 |
| Integrin alpha V beta 3 | Independent sequence partition | 24 | 0.750 | 0.898 |
| Integrin alpha V beta 3 | Matched exact-only exclusion | 25 | 0.600 | 0.789 |
| Integrin alpha V beta 3 | Family exclusion | 25 | 0.600 | 0.744 |

The independent sequence split changes both the test cases and their prevalence. Its AP difference from the family split therefore does not isolate leakage. The matched exact-only diagnostic keeps the test cases fixed. Its ridge AP difference from the family-aware fit is -0.0143 for TG2 (95% interval -0.0663 to +0.1111) and +0.0457 for integrin (-0.0612 to +0.1292). Neither paired interval resolves a directional effect at the primary threshold.

Matched exact-only exclusion changes ridge AP by -0.0143 for TG2 and +0.0457 for integrin, with both paired intervals spanning zero ([paired differences](results/paired_differences.tsv)).

![Average precision under independent sequence splitting, matched exact-only exclusion and family exclusion](figures/03_random_vs_family.png)

Figure 3. AP and 95% family-bootstrap intervals for all three split definitions. Observational references use the same full-count readouts on the matched test cases. The independent sequence partition has a different case mix.

### Generation and novelty

The Markov model uses read-weighted first-order transitions, add-one smoothing and the training length distribution. Generation seeds are 11, 23, 37, 53, 71. Each target-seed run requests 10,000 unique candidates, and all 10 runs attain that budget. The reported total is the sum of within-run uniqueness; cross-run overlap was not quantified. No exact test sequence or positive test family was recovered.

Recovery is defined by an edit-identity edge to an assayed test sequence at the primary threshold. A positive test family contains at least one sequence with an original Positive label. This measures coverage of known sequence neighborhoods; it does not establish binding by an untested candidate. Per-seed outputs and budgets are recorded in [generation statistics](results/generation_statistics.tsv), [budget curves](results/generation_budget_curves.tsv) and [seed variation](results/seed_variance.tsv).

All generation seeds recover zero positive test families at every retained budget for both targets ([budget curves](results/generation_budget_curves.tsv)).

![Positive experimental test-family recovery across fixed Markov generation budgets](figures/05_generation_budget.png)

Figure 5. Positive test-family recovery at 100, 500, 1,000, 5,000, 10,000 unique candidates per run. All five seed curves coincide at zero for both targets.

Median nearest-training edit identity is 0.600 across TG2 seeds and 0.600 across integrin seeds ([novelty distributions](results/novelty_distribution.tsv)).

![Maximum training-sequence edit identity versus Markov log probability per nucleotide](figures/06_novelty_ranking.png)

Figure 6. The first 500 candidates from each seed, or 2,500 points per target, plotted against their nearest identity to the filtered training pool. The vertical axis is Markov log probability per nucleotide, not an experimental activity measurement. Full novelty distributions are summarized in [novelty tables](results/novelty_distribution.tsv).

### Secondary structure

ViennaRNA 2.7.2 folds the assay construct, including constant regions and excluding the poly(A) tether, at 37 C. Sequence-only models are compared with the same models plus MFE, base-pair count, paired fraction and ensemble diversity. The primary ridge AP changes are +0.0000 for TG2 (95% interval -0.0074 to +0.0084) and +0.0031 for integrin (-0.0207 to +0.0259). These fixed fits show no clear AP improvement. They do not test whether structure is biologically irrelevant.

Adding structure changes primary ridge AP by +0.0000 for TG2 and +0.0031 for integrin, without a resolved paired improvement ([ablation comparisons](results/paired_differences.tsv)).

![Average precision for sequence-only and sequence-plus-secondary-structure models](figures/07_structure_ablation.png)

Figure 7. Sequence-only and sequence-plus-structure AP with 95% family-bootstrap intervals on identical experimental test cases. Paired AP differences are recorded separately in [ablation comparisons](results/paired_differences.tsv). Canonical-RNA folding is a proxy for the modified molecules.

### Resources and verification

The largest sampled local footprint is 0.546 GB, against a 13 GB ceiling. The largest measured process-tree RSS is 239.5 MB; minimum sampled free space is 16.167 GB. The footprint includes the environment, project caches, temporary files and Git data. Sampling can miss brief peaks, and initial RSS measurements unavailable under the first sampler remain in the history.

| Scientific stage | Elapsed s | Peak sampled RSS MB | Peak sampled project GB |
| --- | --- | --- | --- |
| splits | 9.510 | 79.6 | 0.482 |
| baselines | 24.918 | 205.0 | 0.481 |
| secondary | 29.962 | 44.0 | 0.481 |
| ablation | 25.520 | 211.4 | 0.481 |
| candidates | 48.516 | 145.9 | 0.481 |
| evaluation | 25.406 | 131.3 | 0.481 |

These timings are the most recent successful executions that updated the scientific stage manifests; cached validation and skip times are excluded. Complete attempts, including failures and restarts, are preserved in [resource records](logs/resource_usage.tsv) and [failure records](logs/failures.tsv).

Across logged attempts, sampled RSS stays at or below 239.5 MB and project size at or below 0.546 GB ([resource records](logs/resource_usage.tsv)).

![Latest per-stage elapsed times and maximum sampled memory and disk across execution attempts](figures/08_resources.png)

Figure 8. Latest elapsed time per stage, which includes cached reruns, alongside maximum sampled RSS and project footprint across attempts. These elapsed times are not cold-start fitting costs. The [full-size PNG](figures/08_resources.png) and [PDF](figures/08_resources.pdf) preserve readable stage labels.

The recorded test suite has 54 passing tests, 0 failures and 0 skipped tests. Clean-clone smoke and test checks, shellcheck, restart validation, numerical traceability and leakage checks are recorded in the [verification audit](results/verification.tsv). Full-mode refusal is documented separately from passing local checks.

## Repository structure

```text
config/       benchmark decisions, primer design, assay and split manifests
scripts/      numbered stages, resource controls, analysis and reporting
remote/       RaptGen job wrapper and AptaDiff execution blocker
results/      compact measurements, predictions, provenance and uncertainty
figures/      PNG and PDF outputs generated from result tables
logs/         resources, failures, restart manifests and verification
data/         ignored sequencing, external, generated and remote artifacts
tests/        synthetic workflow, scientific checks and artifact validation
docs/         methods, score definitions and analysis report
```

## Usage

16 GB RAM, maximum 13.0 GB local project footprint, with a 1 GB external free-space reserve. Run from the repository root with Git, Bash 4 or later and a compatible Python interpreter. The recorded Windows environment used Python 3.13.11, Git Bash 5.2.37 and the versions in [requirements-lock.txt](config/requirements-lock.txt). The local environment contains no GPU PyTorch installation or neural checkpoint.

For a fresh installation, set `BOOTSTRAP_PYTHON` to the intended Python executable if the interpreter on `PATH` is unsuitable. For example, on a Linux or WSL installation with Python 3.13 available:

```bash
BOOTSTRAP_PYTHON=python3.13 bash scripts/02_install.sh
bash scripts/00_configure.sh --threads 1 --ram 16 --disk 13 --gpu-mode none --yes
bash run_all.sh --mode smoke
bash run_all.sh --mode core
bash scripts/18_verify.sh
```

The smoke workflow uses synthetic reads and requires no public-data download or GPU. The core workflow retrieves both public datasets, runs the local analysis and writes the reports. To resume from prepared upstream outputs:

```bash
bash run_all.sh --mode core --from 10
```

`--target tg2` or `--target integrin` limits retrieval; the comparison still requires processed data from both targets. `--from` requires the preceding validated inputs. Primary data acquisition is independent of the synthetic clean-clone check; final verification reuses provider-validated compact counts rather than downloading every FASTQ again.

For hosted RaptGen work, prepare the jobs locally, then copy the selected bundle and tracked remote wrapper to the host:

```bash
# Local preparation; measured stage resources are in pipeline_stage_summary.tsv.
source scripts/lib.sh
stage gpu_jobs scripts/10_prepare_gpu_jobs.py
# On the hosted Linux GPU machine:
conda env create -f remote/raptgen/environment.yml
conda activate raptgen_published
bash remote/raptgen/run_remote.sh remote/bundles/raptgen_tg2_family output/tg2_family
# After copying the listed return files back to the local repository:
bash scripts/import_remote.sh --directory remote/returned/tg2_family \
  --job remote/bundles/raptgen_tg2_family/job.json
```

Repeat for both targets and the matched diagnostic as described in [remote execution](remote/README.md). These hosted commands are prepared and have no measured production run. After returned artifacts, licensing and neural aggregation are resolved, full evaluation would run with:

```bash
bash run_all.sh --mode full --gpu-mode hosted --from 11
bash scripts/18_verify.sh
```

Full mode currently refuses those unresolved requirements. [HANDOVER.md](HANDOVER.md) lists them explicitly. The analysis is also available as [Markdown](docs/analysis_report.md) and [HTML](results/report.html); Quarto is optional.

## Limitations

Only two targets are evaluated. TG2 discrimination is the least certain result: one negative sequence and 6 independent families provide little information about false-positive ranking. The published assay panels are selected candidates, not random samples of the entire sequence space. Bootstrap intervals describe the fixed test panels and fits; they exclude training-set and hyperparameter uncertainty. Family definitions are arbitrary, and pool capping restricts the neighborhood graph. At the 70% sensitivity threshold, integrin's grouped assay-supervised fit has only one training class; that fit and its structural variant are excluded and logged.

Amplification and selection bias can enter both learned scores and generated sequences; generative log probability measures resemblance to that distribution rather than affinity. Canonical-RNA parameters do not reproduce 2'-fluoro chemistry exactly. No generated candidate has a new wet-lab assay. Hosted GPU dependence and legacy neural environments remain unvalidated, and neural evaluation is unfinished. AptaDiff's paper declares MIT while its inspected repository has no licence file; InstructNA likewise lacks clear repository and weight terms. Both are excluded. Tertiary prediction is disabled and supplies no evidence here. These results do not establish reliable aptamer design without target-specific selection or functional data.

## Data availability

Public sequencing accessions are [DRA009383](https://www.ebi.ac.uk/ena/browser/view/PRJDB9110) and [DRA009384](https://www.ebi.ac.uk/ena/browser/view/PRJDB9111), retrieved on 2026-10-03 (Asia/Taipei). DDBJ-submitted INSDC metadata is resolved through the ENA Portal API and experiment XML; FASTQs use the exact provider URLs in [run metadata](results/selex_manifest.tsv). That table records provider checksums, expected sizes, read counts and round-assignment evidence. [Download provenance](results/download_manifest.tsv) records local hashes and timestamps. The primary supplement is retrieved from the publisher and parsed directly. `bash run_all.sh --mode core` regenerates these inputs and derived results.

Compact results and derived assay records are tracked. Raw FASTQs, compact sequencing counts, original PDFs, complete generated candidate tables, environments and model weights are ignored by Git and regenerated through the Bash workflow. Retained result tables support inspection without downloading the sequencing inputs. Resource summaries are snapshots of the log rows available when the report started; [run_summary.tsv](results/run_summary.tsv) records the snapshot size.

## Citation

The dataset and experimental evaluation originate from [RaptRanker](https://doi.org/10.1093/nar/gkaa484). Model references are [RaptGen](https://doi.org/10.1038/s43588-022-00249-6), [AptaDiff](https://doi.org/10.1093/bib/bbae517) and the excluded optional method [InstructNA](https://doi.org/10.1038/s43588-026-00965-3). Upstream commits, source checks and exclusions are recorded in [citations](results/citations.tsv) and [software provenance](results/software_manifest.tsv).

Software references are [ViennaRNA](https://www.tbi.univie.ac.at/RNA/), [NumPy](https://numpy.org/citing-numpy/), [SciPy](https://scipy.org/citing-scipy/), [scikit-learn](https://jmlr.org/papers/v12/pedregosa11a.html), [Matplotlib](https://matplotlib.org/stable/project/citing.html) and [RapidFuzz](https://github.com/rapidfuzz/RapidFuzz). Family construction uses this repository's exact connected-component algorithm and RapidFuzz global edit distances.

## Licence

Repository code is released under the [MIT licence](LICENSE), copyright Qasim Hussain. Dependencies, upstream implementations, datasets and pretrained weights retain their own terms; this licence does not replace them. RaptGen is MIT and is installed separately on the hosted path. ViennaRNA has a custom licence with attribution and redistribution conditions. No pretrained weights are redistributed.

The RaptRanker article and supplement carry CC-BY-NC terms. The original PDF and sequencing inputs are not redistributed here. Derived factual assay fields are attributed to the supplement, and [data terms](results/data_terms.tsv) are recorded separately from the code licence. Software and data terms were checked on 2026-10-03; installed dependency terms are recorded in [licence metadata](results/transitive_licences.tsv).
