# COB-GHIST
GHIST challenge with Chicago Open Bioinformatics

## Installation & Setup

### Python: Retrieving Data Files

Data files are pulled from Synapse programmatically via [pysrc/synapse_retrieval.py](pysrc/synapse_retrieval.py). To use it:

1. Create a Python environment (3.10–3.13; `synapseclient` doesn't yet support 3.14+) and install the project from the repo root:

   ```bash
   conda create -n ghist26 python=3.12
   conda activate ghist26
   pip install -e .
   ```

2. Add a Synapse [personal access token](https://help.synapse.org/docs/Managing-Your-Account.2055405596.html#ManagingYourAccount-PersonalAccessTokens) to a `.env` file in the repo root (gitignored):

   ```
   SYNAPSE-PAT=<your-synapse-personal-access-token>
   ```

3. Run the script to download and unpack the data into `data/`:

   ```bash
   python pysrc/synapse_retrieval.py
   ```

### R (lang)

### Nextflow

After running the survey pipeline, generate an offline plotting report from the
repository root:

```bash
pip install -e '.[plotting]'
python -m pysrc.plot_results
```

Open `nextflow/results/plots.html`. See the [report documentation](nextflow/README.md#plotting-results)
for supported formats, units, and custom results directories.

The report includes a combined graphic with a key for each challenge, demographic
time scaled to generations when published assumptions allow it, and explicitly
exploratory target estimates. Downloadable PNGs are saved under
`nextflow/results/plot_summaries/`; estimates and provisional pair/region
assignments are exported to `nextflow/results/challenge_estimates.json`.
Targets unsupported by the available outputs are marked unavailable.

Python sources are rooted at `pysrc/`: setuptools discovers packages there and
installs the `plot_results` and `synapse_retrieval` modules. With an editable
installation, `python -m plot_results` also runs the report generator.

## Survey of Existing Tools

### SweepFinder2 / Evo 2 experiments

The independent `sweep-experiment` runner supports SweepFinder2 map comparisons,
optional local Evo 2 context-offset scoring, reference-anchored neutral msprime
controls, simulation-based B estimates, and held-out null calibration.
Choose PyTorch binaries with `uv sync --python 3.12 --extra experiments --extra rocm`
(AMD), `--extra cuda` (NVIDIA), or `--extra cpu`. Select only one backend;
omit the backend extra for runs without PyTorch.
See [setup, manifest format, examples, and hardware limitations](docs/sweep_experiments.md).

You may need to install nextflow: [Install Instructions](https://docs.seqera.io/nextflow/install#self-install). This makes it possible to quickly run each of the legacy tools in their own environment rather than lengthy configuration for local installation. Some of these tools are also quite old and may have maintenance issues making them incompatible with modern hardware / ecosystems. Using docker containers, nextflow can isolate and run each tool in its own environment.
