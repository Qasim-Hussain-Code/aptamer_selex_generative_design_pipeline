# AptaDiff execution blocker

The inspected commit is ff2719c7840c18d2b912fd9fa4953ca5bd0f5602. The primary paper's code-availability statement declares MIT, but the complete repository tree contains no licence file and GitHub reports no licence. The separate model-weight and dataset terms are also unresolved. This repository takes the conservative exclusion required by the benchmark specification. No AptaDiff source is vendored, modified or executed.

The data files were inspected as evidence. Dataset A is IGFBP3 and dataset B is PTK7. The primary paper's Table 1 assigns the public TG2 and integrin data to C and D. Its availability sentence associates the DRA accessions with A and B. The contradiction is retained in results/data_provenance.tsv.

After an explicit software licence and weight terms are supplied by the authors, re-run the source inspection, update the recorded commit and licence, and implement the published VAE-to-conditional-diffusion interface in this isolated directory. The current launcher deliberately exits with status 78. Its environment file is a preflight placeholder, not a runnable claim of reproduction.

Resume preflight with `bash scripts/acquire_reference_inputs.sh --part upstream` and `bash run_all.sh --mode core --from 01`. Full mode remains blocked until this specific implementation gap and the missing returned artifacts are resolved.
