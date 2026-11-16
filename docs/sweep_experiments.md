# SweepFinder2 and Evo 2 experiments

This is an exploratory benchmark, not a validated sweep caller. The console
script `sweep-experiment` (also `python -m sweep_experiments`) provides `run`,
`simulate`, `estimate-b`, and `evaluate` subcommands. The modules are independent
of the existing Nextflow survey and do not modify its outputs.

## Environment and reference sequences

Use Python 3.11 or 3.12 for Evo 2, following the upstream runtime requirements.
Install SweepFinder2 and bcftools in your runtime environment (for example via
Bioconda), then:

```bash
uv sync --python 3.12 --extra experiments --extra rocm
uv run --no-sync sweep-experiment --help
```

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) first.
Select **one** PyTorch backend:

| Backend | Installation | PyTorch repository |
|---|---|---|
| AMD / Linux | `uv sync --python 3.12 --extra experiments --extra rocm` | ROCm 6.3 |
| NVIDIA / Linux or Windows | `uv sync --python 3.12 --extra experiments --extra cuda` | CUDA 12.8 |
| CPU / Linux, Windows, macOS | `uv sync --python 3.12 --extra experiments --extra cpu` | CPU index; PyPI native wheel on macOS |
| No PyTorch (SF2/simulation only) | `uv sync --python 3.12 --extra experiments` | None |

These extras pin PyTorch 2.7.1, an upstream-recommended Evo 2 baseline, rather
than silently taking a newer binary/kernel ABI. The backend version names refer
to the wheel runtimes; confirm your GPU and installed driver are compatible.
CPU PyTorch does not enable CPU inference in the Evo 2 adapter.
The package manager does not auto-detect hardware. Unsupported platform/wheel
combinations fail to resolve rather than silently selecting a different backend.

The extras are mutually exclusive under uv; do not use `--all-extras`.
Switch backends by repeating `uv sync` with only the desired backend extra.
The checked-in `uv.lock` records all backend branches; add `--locked` to require
that exact resolution without updating it.
`uv sync` manages the project's `.venv`, not an activated Bioconda environment:
keep external bcftools/SweepFinder2 on PATH or specify their executable paths.
Use `uv run --no-sync ...` after syncing, or pass the same extras to `uv run`,
so a default sync does not remove the selected optional dependencies.

Only explicitly routed packages use the PyTorch indexes; ordinary dependencies
still come from PyPI. Plain pip ignores `[tool.uv.sources]` and uv's extra
conflicts. For a manually managed pip/Conda Python environment instead, install
PyTorch from exactly one official index, then the project **without** a backend extra:

```bash
# AMD example; substitute /cu128 or /cpu for another backend.
python -m pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/rocm6.3
python -m pip install -e '.[experiments]'
```

The local `rustysweepfinder` directory is currently a standalone Cargo crate,
not an installable Python package. The former unresolved `rssf` Python extra
has been removed; build the Rust crate with Cargo separately until bindings
and Python build metadata exist.

Verify the installed binary in the selected environment:

```bash
uv run --no-sync python -c "import torch; print(torch.__version__, torch.version.cuda, torch.version.hip); print('GPU available:', torch.cuda.is_available())"
```

The optional extra supplies `biocommons.seqrepo`, msprime, and Biopython.
SweepFinder2 and bcftools remain external executables. Evo 2 is intentionally
not installed by this extra: first install a hardware-compatible PyTorch and
the upstream Evo 2/Vortex kernel dependencies, then the official `evo2` package
using that environment's Python. Auxiliary installs are not managed by the
project's extras: keep them in a dedicated environment or reapply them after an
exact `uv sync` (which removes undeclared packages). Do not install an incompatible
PyTorch build as a side effect of installing kernels.
Do not assume that installing a ROCm PyTorch wheel also ports Vortex, FlashAttention,
or Transformer Engine. No sequence data are sent to a remote inference service.
The official Evo 2 loader may download model weights to the Hugging Face cache.

Populate a local SeqRepo with the primary chromosomes of **GRCh38.p14**, the
current GRCh38 patch release checked when this draft was written. Reference
lookup uses versioned accessions, not mutable chromosome aliases:

| Manifest dataset | Chromosome | RefSeq accession | Length |
|---|---|---|---|
| `testing` | 7 | `NC_000007.14` | 159,345,973 |
| `final` | 9 | `NC_000009.12` | 138,394,717 |

Point `--seqrepo` at the local SeqRepo root. Reference sequence provisioning is
explicit; the runner never downloads an entire human reference automatically.
These accessions describe primary chromosomes, not patch scaffolds.
Every scored SNP's REF is checked against SeqRepo before loading Evo 2.
**The GHIST VCFs are simulations: matching chromosome names and lengths do not
prove that their nucleotide sequences match GRCh38. A mismatch is an error,
not an invitation to substitute reference bases.** SF2-only runs do not need SeqRepo.

## Run existing VCFs

```bash
sweep-experiment run \
  --manifest examples/sweep_experiments.json \
  --output experiment_results/sf2-baseline \
  --configurations physical,constant --grid 1000
```

`--sweepfinder` and `--bcftools` override executable names or paths.
Output directories must be new. Failures propagate as exceptions with command
logs retained; a run is finished only when `run.json` says `complete`.

The example manifest intentionally labels the challenge `unknown` and
`exploratory`: it is not independent simulation truth. Paths inside manifests
are relative to the manifest, not the current working directory. Each case
requires `id`, `vcf`, `mask` (callability BED), `dataset`, `scenario`, and `split`.
Use one sorted contig per case. Scenario values are `neutral`, `bgs` (no sweep),
`sweep`, `bgs_sweep`, or `unknown`; splits are `calibration`, `evaluation`, or
`exploratory`. Set `sweep_position` to a 1-based position when known.

The bcftools adapter streams genotype queries, excludes duplicate-position
records, multiallelic/non-SNP records and sample-monomorphic sites, and counts
called allele copies. It supports missing calls and does not assume REF is
ancestral. Folded mode is the default; `--unfolded` requires valid `INFO/AA`.
Masks use 0-based half-open BED intervals. VCF and SF2 positions are 1-based.
No variant quality filters beyond these are invented; perform your desired QC
before constructing the manifest.

Each case produces:

- Per-site `.freq` input and the actual `SweepFinder2 -f` background `.spect`.
- SF2 scans and individual subprocess logs.
- Optional Evo per-SNP/per-context/per-strand probability rows.
- Top-level `summary.tsv`: maximum score, peak coordinate, truth, localization error.
- Top-level `run.json`: arguments, input SHA-256 digests, commands, reference
  accessions and selected model/backend metadata when Evo is enabled.

By default the background SFS is estimated from the case itself. This is the
self-background baseline, not an independent neutral reference. A case may
instead specify either:

- `background_vcf` plus `background_mask`: estimate an SFS from separate sites.
- `background_spectrum`: an externally prepared, compatible SF2 spectrum.

Calibration and evaluation VCF paths must be disjoint; no evaluation VCF may be
used as an external background. This cannot detect copied datasets or shared
ancestry between separately named files: independent replicates and appropriate
background spectra remain the experimenter's responsibility.

## Recombination and B maps

Configurations:

| Name | Mode | Additional case fields |
|---|---|---|
| `physical` | SF2 `-l`, no map | None |
| `constant` | SF2 `-lr`, constant rate | None; `--constant-rate` is per base per generation |
| `map` | SF2 `-lr`, spatial map | `recombination_map` |
| `bmap` | SF2 `-lrb`, spatial map and B correction | Fields below |

Supply a cumulative genetic map TSV with exactly named columns:

```text
position	cM
1	0
100001	0.1
200001	0.3
```

It must span every SF2 input position and have strictly increasing physical
positions and nondecreasing finite cM. Linear interpolation followed by
differencing produces SF2's **interval distances in cM**, first row zero.
This is not a raw rate table or a cumulative map passed directly to SF2.
Distances span inaccessible gaps; masking does not collapse coordinates.
The map is a single-contig map, interpreted for that case's VCF contig.

For `bmap`, also supply:

```json
{
  "sf2_frequency": "prepared-polymorphism-and-fixed-differences.freq",
  "background_spectrum": "independent-compatible-background.spect",
  "recombination_map": "recombination.tsv",
  "bmap": "background.bedgraph",
  "n1": 10000,
  "n2": 10000,
  "divergence_time": 100000
}
```

Run with `--unfolded`. `sf2_frequency` has the normal
`position x n folded` tab-delimited header; it must include matching counts for
all eligible VCF SNPs and genuinely observed fixed derived differences (`x=n`).
Every position must be callable. External frequency input is used for **all**
SF2 configurations in that run, so configuration comparisons use the same sites.
It is the user's responsibility to prepare reliable outgroup alignment and a
background spectrum with matching ascertainment/sample-size support.
The runner does not infer substitutions from SNP absence or fabricate outgroup data.

`n1`, `n2` and divergence time are current ingroup effective size, ancestral
effective size, and ingroup/outgroup split time in generations. The correction
assumes a suitably diverged outgroup; this runner requires `T/(2*N1) >= 4`,
following the manual's approximate guideline. This ratio alone is not proof of
ingroup monophyly or a good demographic model.

B maps are headerless `chrom start end B` bedGraph: sorted, nonoverlapping,
0-based half-open intervals covering every SF2 input position, with `0 < B <= 1`.
Use Bvalcalc or another justified independent model to generate them.
Both B and recombination maps are sampled at **frequency-input positions**,
not scan-grid positions.

## Evo 2 context-position experiments

```bash
sweep-experiment run \
  --manifest examples/sweep_experiments.json \
  --output experiment_results/evo-contexts \
  --configurations physical,constant \
  --evo --seqrepo /path/to/seqrepo \
  --model auto --allow-rocm \
  --context 8192 --offsets 1024,4096,7168 --evo-stride 1000
```

`auto` considers 40B, 7B, then 1B, skipping checkpoints incompatible with
documented hardware/context requirements or conservative available-memory
estimates. Each attempted load is probed with the requested context length.
Only out-of-memory errors cause a smaller-model retry; missing kernels, broken
imports, and numerical failures propagate rather than being hidden as a
successful downgrade. A later real-sequence inference OOM is also an error;
the runner does not mix checkpoints within a result set.

Upstream currently documents NVIDIA FP8/Transformer Engine requirements for
40B and 1B, but BF16 support for 7B. Consequently, on an AMD 7900 XTX this
adapter can attempt **7B only**, with explicit `--allow-rocm`; actual support
depends on the installed Vortex/kernel stack and is not guaranteed. It does
not implement a pure-PyTorch port, quantization, or CPU offload. Multi-GPU
placement is left to upstream Vortex. GPU fit estimates are screening heuristics,
not guarantees of maximum possible fit; the context probe tests actual allocation.

For each sampled SNP, the window slides so the target appears at each requested
0-based offset, on both forward and reverse-complement strands. The code reads
logits at `target - 1` to obtain causal next-token probabilities at the target.
Offsets must be >=1; unavailable edge contexts are omitted and wholly
unscorable sites get an explicit `contig_edge` row. Contexts containing ambiguous
bases get `ambiguous_context` rows.

Because Evo 2 is causal, tokens **after** the target cannot influence its
prediction. These experiments vary the amount of preceding context; reversing
the strand probes the other genomic flank. These are reference contexts,
not reconstructed sample haplotypes, and do not infer sweep age or strength.

The output contains log probabilities of REF and ALT, their log-probability
difference, and A/C/G/T probabilities in genomic orientation. Probabilities are
from the full model vocabulary, not renormalized over four bases. Scores use
natural logarithms. Reverse-complement probability columns are complemented
back to genomic orientation.

As an intentionally simple comparator, summary scores average `abs(delta_logp)`
across contexts/strands per SNP, then average per fixed physical window and take
the maximum window score. `--window-size` controls those bins. Sparse bins can
be noisy; this is a functional-compatibility proxy, **not a population-genetic
likelihood or a calibrated sweep probability**. Fix scoring direction,
stride, offsets, and window size before held-out evaluation.

## Neutral controls and external SLiM datasets

```bash
sweep-experiment simulate \
  --seqrepo /path/to/seqrepo --dataset testing \
  --start 10000000 --end 11000000 \
  --replicates 100 --samples 25 \
  --ne 10000 --mutation-rate 1.2e-8 --constant-rate 1e-8 \
  --seed 1729 --output experiment_results/neutral-controls

sweep-experiment run \
  --manifest experiment_results/neutral-controls/manifest.json \
  --output experiment_results/neutral-scans \
  --configurations physical,constant,map
```

Choose an ACGT-only interval. `--dataset final` switches to chr9.
`--recombination-map` can replace the constant rate with a supplied cumulative
map. The simulator uses msprime diploid ancestry and JC69 mutations, conditioning
each site's ancestral nucleotide on the reference through a nucleotide
permutation. JC69 symmetry makes this valid for this simple mutation model.
VCFs include true ancestral alleles and actual GRCh38 REF bases.
This is **not** a context-dependent human mutation model or a guarantee of
realistic language-model score distributions. Start with small intervals.

The output includes VCFs, a callability mask, a cumulative genetic map, and a
manifest. The first half of independent replicates is calibration, the second
half evaluation. Simulation parameters and msprime version are recorded.
Neutral-only controls cannot measure sweep sensitivity or BGS robustness.

Generate BGS-only, sweep-only, and BGS+sweep datasets separately with SLiM and add
them as cases. Record known sweep positions; retain no-sweep BGS controls.
Match sample size, demographic and mutation models, region length, masking,
and scoring configuration across comparisons. A BGS simulation is not reproduced
merely by lowering a constant msprime population size.

### Estimate B from matched simulations

Compute diversity at the **same neutral callable sites** in independent neutral
and BGS-only replicates, not at deleterious sites or in sweep simulations.
Supply a TSV:

```text
scenario	replicate	chrom	start	end	pi
neutral	1	7	0	100000	0.001
bgs	1	7	0	100000	0.0005
neutral	2	7	0	100000	0.0012
bgs	2	7	0	100000	0.0007
```

```bash
sweep-experiment estimate-b --diversity diversity.tsv \
  --bootstrap 1000 --seed 1729 --output b-estimates.tsv
```

The estimate is `mean(pi_BGS) / mean(pi_neutral)`, with paired replicate
bootstrap 95% percentile intervals. Matching IDs denote matched parameter
replicates; at least two are required. No zero-diversity denominator is silently
replaced. Estimates above 1 are reported, not clipped: they may reflect sampling
noise or a misspecified comparison. Investigate them before constructing a
valid SF2 B map. Select `chrom/start/end/B` columns (without a header) to create a
bedGraph only after validating the estimates; retain uncertainty separately.

## Held-out evaluation

```bash
sweep-experiment evaluate \
  --summary experiment_results/neutral-scans/summary.tsv \
  --alpha 0.05 --output experiment_results/evaluation.json
```

Calibration uses only `neutral`/`bgs` no-sweep cases, separately for each score
configuration. Each case contributes its maximum score. Held-out p-values are
`(1 + count(null_max >= observed_max)) / (1 + number_of_null_replicates)`.
Ties are conservative. At least 19 calibration cases are needed to resolve
alpha=0.05, and substantially more are recommended for stable estimates.
The JSON reports per-case detections, localization errors, and detection rates
by truth scenario: rates for neutral/BGS are false-positive rates, while
sweep/BGS+sweep rates estimate power.

Null strata are pooled as supplied: use separate manifests/runs for different
demographies, region lengths, ascertainment schemes, or other incompatible nulls.
This controls a per-region maximum, not the false discovery rate over an
arbitrary number of chromosomes. Do not select score directions, context
positions, windows or maps using evaluation labels.

## Validation and sources

```bash
SWEEPFINDER2_TEST_BINARY=/path/to/SweepFinder2 \
  python -m unittest discover -s tests -p test_sweep_experiments.py -v
```

Tests use real bcftools and, when available, real SF2, SeqRepo and msprime.
GPU-free tests check causal logit indexing, context offsets and strand mapping;
they do not establish numerical equivalence or hardware compatibility of Evo 2.
Run upstream Evo 2 numerical tests in the target GPU environment before trusting
ROCm results.

- [Official SF2 manual](https://github.com/mdegiorgio/SweepFinder2/blob/main/Manual_SweepFinder2_v1.0.pdf),
  sections 3.2, 3.3 and 5.3-5.5.
- [Evo 2 runtime and hardware requirements](https://github.com/ArcInstitute/evo2).
- [Evo 2 7B checkpoint](https://huggingface.co/arcinstitute/evo2_7b).
- [GRCh38.p14 primary chromosome accessions](https://www.ncbi.nlm.nih.gov/grc/human/data).
- [Bvalcalc parameters](https://johrilab.github.io/Bvalcalc/guides/params.html).
