# COB-GHIST

GHIST challenge with Chicago Open Bioinformatics.

## Installation & Setup

With [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
installed, run from the repository root:

1. **Create one environment** for Python, SweepFinder2, and bcftools:

   ```bash
   conda create -y -n ghist26 --override-channels -c conda-forge -c bioconda \
     python uv sweepfinder2 bcftools
   conda activate ghist26
   ```

2. **Install the project and PyTorch into it:**

   ```bash
   uv pip install --python "$CONDA_PREFIX/bin/python" -e .
   TORCH_BACKEND=auto  # AMD override: use a compatible rocm channel if needed.
   uv pip install --python "$CONDA_PREFIX/bin/python" --upgrade torch --torch-backend "$TORCH_BACKEND"
   ```

Keep `ghist26` activated and run commands directly: no separate virtual
environment or `uv sync` is needed. All project dependencies, including msprime,
SeqRepo, and Biopython, install by default. `uv pip` does not use `uv.lock`.

uv selects the PyTorch backend automatically. For AMD, override `TORCH_BACKEND` if
needed with a compatible ROCm channel listed by `uv pip install --help`.
GPU drivers must already be installed; detection may fall back to CPU.
**NTv3 scoring requires a working CUDA/ROCm GPU and access to gated model weights.**
It uses the official PyTorch checkpoint, not Evo2/Vortex or Transformer Engine.
See the [NTv3 and reference setup guide](docs/sweep_experiments.md#environment-and-reference-sequences)
for model access, hardware validation, and licence requirements.

## Download data

Copy [example.env](example.env) to [`.env`](.env), then replace the placeholder with your Synapse
[personal access token](https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens):

```bash
cp -n example.env .env
```

Then download and unpack data into [data/](data/). The token file is gitignored.

```bash
python -m synapse_retrieval
```

## Survey pipeline and reports

Install [Nextflow](https://docs.seqera.io/nextflow/install#self-install) and follow
the [pipeline instructions](nextflow/README.md). After the pipeline finishes,
generate the offline report:

```bash
python -m plot_results
```

Open [plots.html](nextflow/results/plots.html). The report includes challenge
summaries, downloadable PNGs, and exploratory target estimates; unsupported
targets are marked unavailable. See the
[report documentation](nextflow/README.md#plotting-results) for output formats
and custom result directories.

## SweepFinder2 / genomic language model experiments

The independent `sweep-experiment` runner supports SweepFinder2 map comparisons,
optional local NTv3 masked-base or HyenaDNA causal next-base scoring, reference-anchored neutral msprime
controls, simulation-based B estimates, and held-out null calibration.
Runs also export GHIST-format **headerless BED3** sweep regions, with empirical
SF2 p-values and genomic language model evidence in separate decision tables. Significant calls
require matched no-sweep calibration controls; otherwise only candidates are
reported. See [BED exports and calibration](docs/sweep_experiments.md#ghist-bed-exports-and-decision-evidence).
Check available commands:

```bash
sweep-experiment --help
```

```
sweep-experiment run \
  --manifest examples/sweep_experiments.json \
  --output experiment_results/sf2-ntv3 \
  --configurations physical,constant \
  --ntv3 \
  --seqrepo /usr/local/share/seqrepo/2024-12-20
```

Use `--hyenadna` instead of `--ntv3` for the public HyenaDNA large checkpoint.
Start in a new output directory with `--glm-stride 1000`; see the
[HyenaDNA example](docs/sweep_experiments.md#hyenadna-causal-context-position-experiments)
for ROCm and window-size guidance. See the [experiment guide](docs/sweep_experiments.md) for manifests, examples,
reference data, and hardware limitations.
