# Data storage

The raw FASTQ files, external PDFs, upstream inspection files, compact sequencing counts, generated sequences and returned GPU artifacts are ignored by Git. All count tables use transcribed-strand variable regions with T encoding. U conversion occurs only for canonical RNA folding. Molecular reverse-complement equivalence is disabled; reverse-oriented sequencing reads can still be oriented by their constant regions.

Reconstruct the external metadata and primary supplement with the Bash commands in README.md. results/download_manifest.tsv records source URLs, retrieval times, provider MD5 values and locally calculated hashes. results/selex_manifest.tsv resolves rounds from experiment library names. scripts/04_fetch_selex.sh verifies counts and removes each raw file before retrieving the next one. No round was read-subsampled in the executed run.

Derived assay measurements are attributed to the primary RaptRanker supplement. The article/supplement terms are recorded separately from this repository's MIT code. No upstream supplementary PDF, model source or checkpoint is redistributed here. Compact result tables preserve the facts needed to review the analysis. Training-pool capping is separate from read processing and is defined in config/benchmark.yml.
