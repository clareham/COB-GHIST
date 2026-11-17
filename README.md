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
**Evo 2 requires separate model/kernel setup and a supported GPU**; installing
ROCm PyTorch alone is insufficient. See the
[Evo 2 and reference setup guide](docs/sweep_experiments.md#environment-and-reference-sequences).

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

## SweepFinder2 / Evo 2 experiments

The independent `sweep-experiment` runner supports SweepFinder2 map comparisons,
optional local Evo 2 context-offset scoring, reference-anchored neutral msprime
controls, simulation-based B estimates, and held-out null calibration.
Check available commands:

```bash
sweep-experiment --help
```

```
sweep-experiment run \
  --manifest examples/sweep_experiments.json \
  --output experiment_results/sf2-evo \
  --configurations physical,constant \
  --evo --allow-rocm --model evo2_7b \
  --seqrepo /usr/local/share/seqrepo/2024-12-20
```

See the [experiment guide](docs/sweep_experiments.md) for manifests, examples,
reference data, and hardware limitations.
