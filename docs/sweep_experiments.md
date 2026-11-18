# SweepFinder2 and NTv3 experiments

This is an exploratory benchmark, not a validated sweep caller. The console
script `sweep-experiment` (also `python -m sweep_experiments`) provides `run`,
`simulate`, `estimate-b`, `evaluate`, and `export-regions` subcommands. The modules are independent
of the existing Nextflow survey and do not modify its outputs.

## Environment and reference sequences

Follow the [two-step environment setup](../README.md#installation--setup).
All experiment dependencies are installed by default; there are no experiment
or backend extras to select. The README's `TORCH_BACKEND` setting is the only
backend switch needed.
The project does not pin PyTorch versions or CUDA/ROCm repository channels.
Both the project and PyTorch are installed by `uv pip` into the activated
Conda environment; no additional virtual environment is created.
`--torch-backend auto` detects installed GPU drivers and chooses a compatible
official PyTorch index. If no supported GPU is detected, it falls back to CPU;
check the selected build and actual GPU execution below before loading NTv3.

To select CPU explicitly, replace `auto` with `cpu`. To override GPU detection,
replace `auto` with a CUDA or ROCm backend identifier supported by your installed
uv release (see `uv pip install --help`). Runtime-channel selection then lives
in your installation command, not in this project's configuration. Upgrade uv
if it does not recognize a newly published backend; do not assume detection
supports every GPU or distribution. If uv still lacks that backend, install
from the official index selected by the PyTorch installer instead:

```bash
uv pip install --python "$CONDA_PREFIX/bin/python" --upgrade torch \
  --index-url <official-index-url>
```

When changing backends in an existing
environment, add `--reinstall-package torch` to ensure the build is replaced.
PyTorch selects its matching kernel dependencies, including ROCm Triton when
required; the project does not independently constrain Triton.

Different indexes may offer different PyTorch versions, so do not attribute
score differences solely to hardware. Automatic selection does not validate
NTv3 inference on your GPU. The adapter uses the official NTv3 PyTorch model;
no Evo2, Vortex, Transformer Engine, or separately installed FlashAttention
package is required by this adapter. Run numerical checks on the target
hardware before interpreting results.

The backend version names refer to the **wheel build runtimes**, not the installed
ROCM-SMI library. For example, `ROCM-SMI-LIB version: 7.8.0` identifies that
management component; it neither proves the HIP runtime version nor implies
that a `rocm7.8` PyTorch wheel repository exists.
Confirm your GPU, kernel driver, and host OS
are compatible with the selected build; a newer management-library version alone
does not establish compatibility.
CPU PyTorch does not enable CPU inference in the NTv3 adapter.
The documented `uv pip` installation does not consume `uv.lock`. Record
`uv pip freeze --python "$CONDA_PREFIX/bin/python"`, `conda list --explicit`,
and the selected backend with experiment results for reproducibility.
Repeat the PyTorch installation command with `--upgrade` when deliberately
refreshing the backend. Ordinary project dependencies come from PyPI.
Keep the Conda environment activated so Python, SweepFinder2, and bcftools
remain on PATH, and run commands directly without `uv run` or `uv sync`.
The [official PyTorch installer](https://pytorch.org/get-started/locally/)
can help select a wheel repository compatible with your hardware.

The local `rustysweepfinder` directory is currently a standalone Cargo crate,
not an installable Python package. The former unresolved `rssf` Python extra
has been removed; build the Rust crate with Cargo separately until bindings
and Python build metadata exist.

Verify the installed binary in the selected environment:

```bash
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.version.hip); print('GPU available:', torch.cuda.is_available())"
```

For AMD expect a ROCm wheel build, a nonempty `torch.version.hip`,
and `GPU available: True`. To test actual execution, not just enumeration:

```bash
python -c "import torch; assert torch.version.hip; x = torch.ones((256, 256), device='cuda'); print((x @ x)[0, 0].item())"
```

PyTorch intentionally uses the `cuda` device/API names on ROCm too. The matrix
operation should print `256.0`; successful execution still does not validate
the NTv3 model.

The default dependencies include `biocommons.seqrepo`, msprime, and Biopython.
SweepFinder2 and bcftools remain external executables. Transformers is a project
dependency; PyTorch is installed separately using the hardware-compatible command
above. No additional NTv3 Python package or TE source build is needed.

Request access to [InstaDeepAI/NTv3_pretrained](https://huggingface.co/InstaDeepAI/NTv3_pretrained)
and accept the **InstaDeep Open Model Licence (non-commercial use, including
outputs)**. Approval is manual. Authenticate in the same Conda environment:

```bash
hf auth login
```

The loader downloads weights, tokenizer, and custom Python modelling code to the
Hugging Face cache. It uses `trust_remote_code=True`, pinned by default to Hub
commit `fd58207244d828957e633c591f9babd8bb2b970e`; review upstream code before use.
`--model-revision` overrides that revision for all three artifacts. Record and
prefer a commit hash over a moving branch. Access errors and missing dependencies
are reported directly, not replaced with another model.
No sequence data are sent to a remote inference service.
Existing Evo2/Transformer Engine installations are not used by this runner.

Populate a local SeqRepo with the primary chromosomes of **GRCh38.p14**, the
current GRCh38 patch release checked when this draft was written. The runner
explicitly resolves each chromosome's versioned RefSeq accession from
`GRCh38.p14` assembly metadata, then queries that accession in SeqRepo's `NCBI`
namespace with exact 0-based, half-open `start`/`end` coordinates:

| Manifest dataset | Chromosome | RefSeq accession | Length |
|---|---|---|---|
| `testing` | 7 | `NC_000007.14` | 159,345,973 |
| `final` | 9 | `NC_000009.12` | 138,394,717 |

Point `--seqrepo` at the local SeqRepo root. Reference sequence provisioning is
explicit; the runner never downloads an entire human reference automatically.
The repository must contain the corresponding `NCBI` accessions above. The
`2024-12-20` snapshot contains these sequences but lacks `GRCh38.p14` chromosome
aliases, so the runner does not depend on those aliases or fall back to another
assembly.
These accessions describe primary chromosomes, not patch scaffolds.
Every scored SNP's chromosome and 1-based position bounds are checked before
loading NTv3. Simulated VCF REF/ALT labels are retained and are not required to
match the GRCh38 target base. Fetched contexts are not rewritten to match those
labels. **These scores compare simulated alleles in human reference contexts;
they are exploratory, not validated biological mutation effects.**
SF2-only runs do not need SeqRepo.

Validation against the `2024-12-20` snapshot covered all 4,128,923 eligible,
masked SNPs in the seven locally available sweep VCFs. Their compressed-file
checksums match the current Synapse version 2 releases, whose version comments
confirm the correction to 1-based positions. Nevertheless, REF matches only
24.88-25.02% of unambiguous GRCh38.p14 bases; shifting positions by either one
base does not resolve the disagreement. The coordinate conversion and 6,342
sampled strand/offset context targets were checked separately and were correct.
REF disagreement is therefore expected and does not block NTv3 scoring.
The target base itself is masked before inference.

## Run existing VCFs

```bash
sweep-experiment run \
  --manifest examples/sweep_experiments.json \
  --output experiment_results/sf2-baseline \
  --configurations physical,constant --grid 1000
```

`--sweepfinder` and `--bcftools` override executable names or paths.
Existing run and per-case directories are accepted. Reruns overwrite generated
files; unrelated files and outputs from configurations not rerun are retained.
This reruns the analysis rather than resuming it. Failures propagate as exceptions with command
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
- Optional `ntv3.tsv` per-SNP/per-context/per-strand probability rows.
- Top-level `summary.tsv`: maximum score, peak coordinate, truth, localization error.
- Top-level `run.json`: arguments, input SHA-256 digests, commands, reference
  accessions and model revision, scoring semantics, precision, and backend
  metadata under `ntv3` when scoring is enabled.
- Automatic `regions/` BED3 predictions and separate SF2/gLM decision evidence;
  see [GHIST BED exports](#ghist-bed-exports-and-decision-evidence).

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

## NTv3 context-position experiments

```bash
sweep-experiment run \
  --manifest examples/sweep_experiments.json \
  --output experiment_results/ntv3-contexts \
  --configurations physical,constant \
  --ntv3 --seqrepo /path/to/seqrepo \
  --context 8192 --offsets 1024,4096,7168 --ntv3-stride 1000
```

The checkpoint is `InstaDeepAI/NTv3_pretrained`, an approximately 703M-parameter
single-base masked language model. Track/annotation checkpoints are not accepted.
Inference uses float32 on GPU 0 (`cuda:0` for both CUDA and ROCm), without
quantization, autocast, CPU offload, or multi-GPU placement. ROCm needs no opt-in
flag. This avoids the Evo2/TE dependency chain, but does **not** establish tested
7900 XTX compatibility.

`--context` must be a multiple of **128**, from 128 through 1,048,576 bases,
because of NTv3's downsampling geometry. Inputs are never silently padded,
cropped, or truncated. The default is 8192; begin with that or a smaller aligned
context on a 24 GB card, adjusting offsets to fit. Long-context support is not
a promise that 1 Mb fits in VRAM. A masked synthetic sequence probes the full
requested context immediately after loading. OOM, import, download, and
numerical errors propagate without switching checkpoints, precision, or device.

For each sampled SNP, the window slides so the target appears at each requested
0-based offset, on both forward and reverse-complement strands. The tokenizer
adds no special tokens, so base and token indices coincide. The target token
is replaced with `<mask>` before the model runs, and logits are read at
**`target`**, not `target - 1`. Offsets may include 0 and must be less than the
context length; duplicates and empty offset lists are rejected.
Unavailable edge contexts are omitted and wholly
unscorable sites get an explicit `contig_edge` row. Contexts containing ambiguous
bases get `ambiguous_context` rows.

NTv3 is **bidirectional**: both flanks of the masked target can influence its
prediction. These experiments vary the allocation of left/right reference
context and test strand sensitivity. These are reference contexts, not
reconstructed sample haplotypes, and do not infer sweep age or strength.

The output contains log probabilities of REF and ALT, their log-probability
difference, and A/C/G/T probabilities in genomic orientation. Probabilities are
from the full model vocabulary, not renormalized over four bases. Scores use
natural logarithms:

```text
delta_logp = log P(ALT | context with target masked)
           - log P(REF | context with target masked)
```

Reverse-complement probability columns are complemented back to genomic
orientation. The target is hidden regardless of whether the VCF REF agrees
with the reference assembly, avoiding target-allele leakage.

### Migrating old Evo2 commands and results

Replace `--evo` with `--ntv3` and `--evo-stride` with `--ntv3-stride`.
Remove `--allow-rocm` and old `--model auto` / `--model evo2_*` arguments;
the default is the NTv3 pretrained checkpoint. Old flags fail explicitly.
Use a **new output directory**: per-site files are now `ntv3.tsv`, provenance
uses `ntv3`, and the summary configuration is `ntv3_mean_abs_delta`.
Historical Evo files are not deleted or relabelled. The TSV columns are
unchanged, but these masked bidirectional probabilities are **not comparable**
to Evo2's causal next-base probabilities. Recompute calibration controls
with the same NTv3 model revision, precision, contexts, and scoring settings;
do not reuse Evo2 null distributions.

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

## GHIST BED exports and decision evidence

Every successful `run` ends by exporting `regions/<case>/` for each case and
SF2 configuration. The [GHIST sweep instructions](https://www.synapse.org/Synapse:syn74803798/wiki/640755),
[submission tutorial](https://www.synapse.org/Synapse:syn74803798/wiki/640752),
and [downloaded example](../data/sweep_submission.example.bed) specify region
predictions as BED. Our submission files match that example exactly:
**three tab-separated columns, no header or extra fields**, 0-based half-open
integer intervals, chromosome `7` or `9` (without `chr`). The updated challenge
VCFs are 1-based; these are explicitly converted, not copied into BED.

For example:

```text
7	18500000	18600000
7	34234000	34454000
```

### Files produced

For each configuration (`physical`, `constant`, `map`, `bmap` as enabled):

| File in `regions/<case>/` | Purpose |
| --- | --- |
| `<configuration>.candidates.bed` | All local SF2 peak support intervals, including uncalibrated candidates. **Not significance-filtered.** |
| `<configuration>.sweeps.bed` | Proposed sweep regions passing the empirical SF2 cutoff. BED3 for GHIST submission review. Empty if none qualify or calibration is unavailable. |
| `<configuration>.sf2.tsv` | All callable scan cells: coordinates, raw CLR, empirical p-value, null count, calibration status, and significance flag. |
| `<configuration>.regions.tsv` | Candidate and predicted-region decisions: coordinates, peak location/CLR/p-value, calibration status, and overlapping raw gLM score summaries. |
| `ntv3.sites.bed` | Unique successfully scored gLM SNPs as one-base BED3 intervals. **Evidence coverage, not sweep predictions.** Empty for SF2-only runs. |

The original `<run>/<case>/ntv3.tsv` remains the **raw gLM output**, unchanged,
including signed REF/ALT log probabilities, `delta_logp`, A/C/G/T probabilities,
strand, offset, context window, and unscored-status rows. It is not replaced
with a ranking or rounded BED score. `regions/metadata.json` records its path
and SHA-256, source scan/mask/summary checksums, null counts, column schemas,
alpha, and the region/p-value definitions. `run.json` links the completed export.
Only treat the export as complete when both `run.json` and
`regions/metadata.json` say `complete`.

The region TSV includes `glm_status`, number of scored SNPs and context/strand
rows, mean absolute delta (average per SNP first, then across SNPs), and maximum
absolute raw delta. `no_scored_sites` and `disabled` are explicit, not zero
evidence. Consult the raw signed scores before making biological decisions.
There is **no automatic gLM threshold, SF2/gLM intersection, or combined p-value**:
NTv3 is a sequence-compatibility proxy, not a trained sweep caller. The default
prediction BED selection uses calibrated SF2 evidence; gLM informs manual review.

### SF2 p-values and null requirements

SweepFinder2 emits a **composite likelihood ratio**, not a calibrated p-value.
We estimate an empirical upper-tail p-value for each CLR against **one maximum
CLR per independent no-sweep calibration replicate**, separately by SF2
configuration:

```text
p = (1 + number of calibration maxima >= observed CLR) / (1 + n_null)
```

This uses conservative ties and never returns zero. Comparing scan cells against
null **region maxima** accounts for scanning within a matched region; it is not
a chi-square approximation or a pointwise test against all correlated null grid
points. Region annotations report the p-value of their maximum CLR, not a
separately fitted interval probability.

By default, controls come from `split: calibration`, `scenario: neutral` or
`bgs` in the same manifest. Calibration cases never get calls or self-calibrated
p-values. Evaluation and exploratory cases may receive calls, without using
their truth labels to choose thresholds. Defaults are `--alpha 0.05` and at
least **19 independent calibration replicates** (prefer many more). With fewer
replicates, available p-estimates are reported as `insufficient_nulls` but no
significant sweep BED rows are emitted. With no controls, p-values are `.`,
status is `uncalibrated`, warnings are logged, and candidate/raw evidence is
still exported. No synthetic p-values or significance claims are invented.

To supply controls from a separate completed scan:

```bash
sweep-experiment run \
  --manifest examples/sweep_experiments.json \
  --output experiment_results/sf2-ntv3-calibrated \
  --configurations physical,constant \
  --ntv3 --seqrepo /path/to/seqrepo \
  --null-summary experiment_results/matched-null-scans/summary.tsv \
  --alpha 0.05
```

`--null-summary` **replaces**, rather than pools with, in-manifest calibration.
Its sibling `run.json` must say `complete`; SF2 binary checksum, grid, folding
mode, and constant recombination rate must match. Only calibration no-sweep rows
are used, never held-out evaluation maxima. Repeated calibration VCF checksums
and shared calibration/prediction VCF content are rejected.

**You must additionally match** demographic/BGS assumptions, scanned region
length, sample sizes, mutation model, accessibility/ascertainment, input frequency
format, background spectrum strategy, and recombination/B-map design. The code
cannot establish biological comparability from these files. A 1 Mb neutral
simulation is not a suitable maximum-CLR null for a full-chromosome scan.
Do not pool different null strata. Testing challenge answers are not null
controls. P-values are not corrected across configurations or cases, and are
not genome-wide FDR estimates. Selecting whichever configuration gives the
smallest p-value requires additional multiple-testing control.

### Region boundaries and re-export

Each SF2 grid location represents a cell extending to the midpoints to its
neighbors (midpoint boundaries rounded down), clipped to the scanned
frequency-input span and intersected with the callable mask. The endpoint cells
stop at the first/last frequency-input base. Fractional grid locations are
mapped to their containing 1-based base; sub-base spacing can create overlapping
one-base cells. Adjacent/overlapping significant cells are merged, without
bridging noncallable gaps. Candidates use local maxima; the leftmost location
represents a flat peak.

These intervals are **grid-resolution support regions**, not confidence
intervals, the swept haplotype's extent, or single-base causal localization.
Sparse grids can give very broad predictions. Calibration tests region-level
presence, not the accuracy of inferred boundaries or each claimed causal locus.
No exon restriction or fixed sweep count is silently inferred; in particular,
the multi-sweep GHIST task's exon-only prior needs explicit scientific review
against the provided annotations.

You can change alpha or add matched null controls **without rerunning SF2 or
NTv3**, after a completed BED-enabled run:

```bash
sweep-experiment export-regions \
  --run-dir experiment_results/sf2-ntv3 \
  --null-summary experiment_results/matched-null-scans/summary.tsv \
  --alpha 0.05
```

This overwrites the enabled configurations' export files and updates export
provenance, leaving raw scans/gLM outputs untouched. Keep the callable mask
unchanged; its checksum is checked. Old runs without case-coordinate metadata
need to be rerun. Reruns do not delete unrelated or disabled configurations'
files: use the configurations and file list in the current export metadata.
Review one chosen configuration's prediction BED per challenge before
submitting; the pipeline never uploads or submits anything automatically.

## Validation and sources

```bash
SWEEPFINDER2_TEST_BINARY=/path/to/SweepFinder2 \
  python -m unittest discover -s tests -p 'test_[ns]*.py' -v
```

Tests use real bcftools and, when available, real SF2, SeqRepo and msprime.
GPU-free NTv3 tests use real CPU tensors with a tiny fixture model to check
target masking, same-position logit indexing, full-vocabulary normalization,
shape/finite-value checks, pinned loading, and error propagation. Runner tests
check context offsets, strand mapping, stride, outputs, and provenance.
They do not load the gated 703M checkpoint or establish ROCm numerical equivalence
or hardware compatibility. Validate the real model on the target GPU before
trusting results.

After model access approval and authentication, opt into a real-checkpoint GPU
smoke test (downloads weights and allocates GPU memory):

```bash
NTV3_GPU_TEST=1 python -m unittest discover -s tests -p test_ntv3.py -v
```

`NTV3_TEST_CONTEXT=128` can reduce the smoke-test window; the default is 8192.
This tests finite masked probabilities and independence from the hidden target
allele, not biological validity or cross-hardware numerical equivalence.

- [Official SF2 manual](https://github.com/mdegiorgio/SweepFinder2/blob/main/Manual_SweepFinder2_v1.0.pdf),
  sections 3.2, 3.3 and 5.3-5.5.
- [NTv3 pretrained checkpoint and PyTorch interface](https://huggingface.co/InstaDeepAI/NTv3_pretrained).
- [NTv3 architecture and releases](https://github.com/instadeepai/nucleotide-transformer/blob/main/docs/nucleotide_transformer_v3.md).
- [NTv3 model licence](https://github.com/instadeepai/nucleotide-transformer/blob/main/MODEL-LICENSE.md).
- [GRCh38.p14 primary chromosome accessions](https://www.ncbi.nlm.nih.gov/grc/human/data).
- [Bvalcalc parameters](https://johrilab.github.io/Bvalcalc/guides/params.html).
