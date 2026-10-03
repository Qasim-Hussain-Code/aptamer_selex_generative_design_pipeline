# Hosted GPU execution

The local machine trains no neural model. RaptGen uses its published architecture, trainer, loss and sampler from a pinned MIT-licensed commit. The wrapper uses length-batched variable regions to preserve indels. It samples the standard normal prior in the harmonized arm. This differs from the paper's GMM-center candidate selection and makes no published-protocol reproduction claim. No activity-guided optimization is performed.

Prepare inputs with `bash run_all.sh --mode core --from 10`. Copy the repository code and one directory from `remote/bundles/` to a Linux GPU host. Bundles contain training sequences, evaluation sequences, family IDs and configuration. They contain no assay labels or responses.

On the host:

```bash
conda env create -f remote/raptgen/environment.yml
conda activate raptgen_published
bash remote/raptgen/run_remote.sh remote/bundles/raptgen_tg2_family output/tg2_family
```

Repeat for both targets and both split strategies. This legacy CUDA stack needs compatible hardware such as a T4. The adapter refuses GPUs beyond the supported generation. Environment resolution and training have not been verified on a hosted machine. These are explicit external tasks, not completed results.

Return only `run_manifest.json`, `heldout_scores.tsv` and the compressed generation tables. The manifest records the checkpoint content hash before deletion. The local importer independently verifies returned files and input/configuration/source hashes; it cannot recompute the deleted checkpoint hash. No weight file is needed for local evaluation.

```bash
bash scripts/import_remote.sh --directory remote/returned/tg2_family \
  --job remote/bundles/raptgen_tg2_family/job.json
bash run_all.sh --mode full --gpu-mode hosted --from 11
```

AptaDiff is excluded pending explicit repository licensing and an authorized adapter. Its launcher stops with an explanation. Local stage 12 currently merges local methods only; neural score joins and candidate aggregation require implementation and validation with returned outputs. Full mode refuses that missing integration explicitly. See HANDOVER.md for the remaining work. No pretrained weights are redistributed. Remote GPU caches and environments belong on the hosted machine and must be reported separately from the measured local footprint.

The length-batch adapter duplicates a singleton length bucket to satisfy upstream BatchNorm; this changes that bucket's effective weight. The score uses the upstream PHMM forward calculation with its finite log-probability floor at the posterior-mean latent point. It is a conditional model quantity, not a latent-integrated likelihood. These harmonization choices have not been assessed in a hosted run.
