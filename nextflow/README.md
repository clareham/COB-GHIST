# GHIST 2026 survey pipeline

Nextflow pipeline for genetic history, sweep detection, and relatedness.
Includes exploratory inference workflows, not scientifically validated results.

## Quick start

Install Nextflow (>=24.04) and Conda, then download the challenge data using the
[repository setup instructions](../README.md#installation--setup).
From the repository root:

```bash
cd nextflow

# Run all enabled challenges and tools
nextflow run . -profile conda

# Resume a previous run (reuses matching successful cached tasks)
nextflow run . -profile conda -resume
```

Defaults: inputs in `../data`, outputs in `results/`. Override paths with
`--data_dir /path/to/data --outdir /path/to/results`, or set
`GHIST_DATA_DIR` / `GHIST_OUTDIR`.

### Conda environments

The `conda` profile resolves every package-list environment from
`conda-forge` and `bioconda` only (`--override-channels`), so whatever a
machine's `~/.condarc` adds - notably `defaults`, which fails
non-interactively under the Anaconda terms of service - cannot break or
alter the solve. Environments are cached in `.conda-envs/` next to the
pipeline (override with `NXF_CONDA_CACHEDIR`) rather than under `work/`, so
clearing the work directory does not rebuild them. Environment creation is
allowed up to 3 hours: Nextflow's default is 20 minutes, which the FlexSweep
(TensorFlow) and mrpast (tsinfer/tsdate/numba) environments can exceed on a
fresh machine - that surfaces as every task of the tool failing with
"Failed to create Conda environment" before any tool code runs.

The [mrpast environment](conf/genetic_history/mrpast.yml) is the one
environment file rather than a package list, and its pins were verified
against the real packages: bio2zarr 0.1.x (needed because tsinfer 0.5.1
requires Zarr 2) breaks on NumPy >=2.4, while Bioconda's Python >=3.11
mrpast builds are compiled against NumPy >=2.4.6. The environment therefore
uses `numpy<2.4` and installs `mrpast` and `tsdate` (which is not on
conda-forge or bioconda) from their PyPI wheels, exactly as mrpast's own
documentation does; every runtime dependency of those two packages is
pinned from conda-forge/bioconda, so pip adds nothing else. No conda-build
step or local channel is needed.

### Limit the run

```bash
# Sweep detection only
nextflow run . -profile conda --run_genetic_history false --run_relatedness false

# Genetic history only
nextflow run . -profile conda --run_sweep_detection false --run_relatedness false

# Inspect resolved configuration without running tasks
nextflow config .
```

To disable a challenge or individual tool, change its `enabled` flag in
`conf/<category>/<challenge>.config`. Missing challenge inputs are skipped
with a warning.

Every (challenge, tool) task is independent and tasks use
`errorStrategy = 'ignore'`: a failing tool is reported and skipped instead
of aborting the whole survey (which, combined with one-task-at-a-time
execution, previously meant the first failure hid every other result).
Check the console summary, `results/execution_report.html`, or
`nextflow log <run name> -f status,name,exit,workdir` for failed tasks.

### VCF normalisation stage

The simulated GHIST VCFs are not "clean" biallelic call sets:
`bottleneck.testing.vcf` carries several distinct simulated sites collapsed
onto one coordinate as separate records (7,678 positions, often with
different REF bases), the other genetic-history VCFs contain triallelic
sites (`ALT=C,G`, genotypes with allele index 2; a handful also occur in the
sweep and relatedness VCFs), and the `##INFO=<ID=AA ...` header line is
unterminated. These broke every tool on `bottleneck`, fastsimcoal2/mrpast on
the other genetic-history challenges, and ADMIXTURE on relatedness2/4
(`plink2 --make-bed` rejects multiallelic variants).

Every VCF challenge therefore first passes through the shared
[NORMALISE_VCF](modules/common/normalise_vcf/main.nf) process
([bin/normalise_vcf.sh](bin/normalise_vcf.sh), bcftools): the header is
repaired, every position that is not exactly one biallelic SNP record is
dropped, and the result is bgzipped/indexed. The dropped positions are
emitted as `<challenge>.excluded_sites.bed` and passed to PSMC, MSMC2 and
fastsimcoal2 so they are removed from callable sequence instead of being
counted as monomorphic (PSMC masks the affected 100 bp window). The BED and
a `<challenge>.normalise.tsv` summary are published under
`results/<category>/<challenge>/normalise_vcf/`; the normalised VCF stays in
the work directory. The `bin/` converters additionally guard against stacked
and multiallelic records themselves, so they also behave on raw VCFs.
`INFO/AA` equals `REF` at every site of the genetic-history VCFs, so the
configured `assume_ref_ancestral` option is exact for this data.

SHAPEIT5 stays disabled because the relatedness VCFs are already phased.
Unsupported catalogue tools also stay disabled; enablement is not automatic.

**Only testing inputs are allowed.** Held-out `*.final.*` scoring files
must never be consumed; [GhistUtils](lib/GhistUtils.groovy) enforces the
testing-file naming guard.

## Configuration

| Location | Purpose |
|---|---|
| [nextflow.config](nextflow.config) | Paths, category switches, execution profiles |
| [conf/base.config](conf/base.config) | Shared process resources and reporting |
| `conf/<category>/<challenge>.config` | Challenge inputs, tool enablement, Conda packages, arguments |
| `workflows/<category>.nf` | Routes enabled challenge/tool pairs |
| `modules/common/normalise_vcf/main.nf` | Shared biallelic-SNP normalisation of every VCF challenge |
| `modules/<category>/<tool>/main.nf` | Shared process implementation per tool |

Tool status is descriptive; `enabled` controls execution:

| Status | Meaning |
|---|---|
| `runnable` | Packaged tool with a documented/standard CLI invocation |
| `runnable-unverified` | Packaged tool; exact invocation still needs verification |
| `runnable-with-caveat` | Runs, but input suitability is limited |
| `survey` | End-to-end exploratory workflow with explicit modeling assumptions |
| `excluded` | Not appropriate for the challenge |
| `stub` | Integration and/or installation required; disabled by default |

PLINK uses 22 diploid autosome pairs for human challenges and 38 for dog
(`autosomes` in each relatedness config). ADMIXTURE's `k: 3` remains a
survey choice, not a validated population count.

Tasks request 2 CPUs, 4 GB, and 2 hours by default; fastsimcoal2, mrpast
and ngsRelate request 8 hours (ngsRelate also 4 CPUs) because their
all-pairs/likelihood work on the larger challenges exceeds 2 hours.
`--max_cpus`,
`--max_memory`, and `--max_time` cap **per-task** requests via
`resourceLimits` (requires Nextflow >=24.04). The local executor also sets
`queueSize = 1`: **only one task runs at a time across all tools and
challenges**, while each tool can still use multiple CPU threads. Task order
is not fixed. This prevents concurrent inference tasks from competing for
your 32 GB RAM; it does not guarantee that an individual task fits in RAM.
Local memory requests are scheduling estimates, not enforced RAM limits.
Keep RAM available for Nextflow, the OS, and other applications.

### Genetic-history assumptions and outputs

Tool-specific `options` live beside `args` in each challenge config.
Simulated variant-only VCFs are explicitly assumed callable across their
header-defined contig lengths; missing genotypes and invalid sites are masked.
For other VCFs, supply a callable BED instead of assuming invariant sites.
BAM preparation uses quality/depth masks and joint calls for MSMC2.

- **PSMC:** one individual's curve, preserving invariant/trailing windows;
  publishes the raw curve and an assumptions JSON, not a fitted split time.
- **MSMC2:** per-contig input, phase ambiguity, separate within-population
  and cross-population fits, and a time-aligned relative cross-coalescence
  curve (`rccr.tsv`). Configured sample prefixes identify populations.
  Times use survey mutation rates of `1.25e-8` (human) or `3e-9`
  (Drosophila) per base/generation; these are assumptions, not challenge truth.
- **fastsimcoal2:** callable-length-aware folded 1D/joint 2D spectra, matching
  haplotype sample sizes, configured size-change or split/migration models,
  and best-likelihood results under `<challenge>.fastsimcoal2_out`.
  Bottleneck/growth currently fit one instantaneous size change, not a
  finite-duration bottleneck or continuous growth trajectory. The observed
  spectrum type is selected with `-m` (folded `_MAFpop0.obs`) or `-d`
  (unfolded `_jointDAFpop1_0.obs`); fsc27 also rejects parameter names that
  are substrings of one another, so the templates use distinct names
  (`NDIP`/`HAPSIZE`, `NPOPA`/`HAPA`, ...).
- **mrpast:** VCF-to-Zarr conversion, per-contig tsinfer/tsdate dated ARGs,
  population maps, and demographic process/solve results under
  `<challenge>.mrpast_out`. Epochs, migration, and rates are survey assumptions.

### FlexSweep training

FlexSweep trains a classifier **for each enabled sweep challenge**, then
extracts empirical features and predicts scan windows. Configurable training
defaults use 1,000 neutral and 1,000 mixed hard/soft sweep simulations,
constant effective size 10,000, and reference-as-ancestral polarization.
No background selection or fitted demographic history is modeled: especially
for BGS challenges, predictions are exploratory and not calibrated discoveries.
The staged accessible BED filters SNPs and low-coverage prediction windows.
Outputs include the model, features, predictions, BED calls, and assumptions.
Training uses the upstream maximum of 1,000 epochs with validation-accuracy
early stopping (patience 5), not a fixed short training run.

CPU training is expensive. The module requests 4 CPUs, 16 GB, and 24 hours,
but the default `max_time` caps it at 8 hours. To allow its full request:

```bash
nextflow run . -profile conda --max_time '24h' -resume
```

### AMD GPU acceleration

FlexSweep remains CPU-only; no ROCm profile is enabled. AMD provides
[TensorFlow 2.18.1 ROCm builds](https://rocm.docs.amd.com/projects/install-on-linux/en/latest/install/3rd-party/tensorflow-install.html),
but CachyOS/Arch is outside the listed
[supported operating systems](https://rocm.docs.amd.com/en/docs-6.4.2/compatibility/compatibility-matrix.html).
Using ROCm here would require replacing the TensorFlow dependency in the
FlexSweep environment (not installing two TensorFlow builds together),
checking the exact GPU/driver compatibility, and validating the full package
stack. AMD's prebuilt container is a possible starting point, not a drop-in
replacement for this Conda workflow. Optional GPU support is deferred rather
than adding an untested custom build. Simulations/features remain CPU work.

### SweepFinder2 resolution

All four sweep challenges use `args: '-l 1000'`: **1,000 scan positions
across the contig**, not 1,000-base spacing. Change the count in the
challenge config to adjust resolution; this is a coarse survey default.
Bare `-l` is invalid: it omits the required grid count and causes exit 255.
`-lg g` specifies spacing instead.
The accessible BED is applied, and spectra are folded by default because the
testing VCFs lack ancestral annotations. Unfolded mode requires valid INFO/AA;
REF is not silently treated as ancestral.

After changing arguments, use `-resume`. SweepFinder2 reruns conversion,
spectrum estimation, and scanning together because they share one process.

## Known limitations

- **Scientific accuracy is unvalidated.** Check methods against simulations
  with known parameters. Relatedness testing answers are publicly available
  through the [GHIST wiki](https://www.synapse.org/Synapse:syn74803798/wiki/640740);
  final scoring data remains held out.
- **Low-coverage inputs:** BAM archives require Synapse challenge-team
  access. The PSMC/MSMC2 BAM workflows are unverified against real reads.
- **Runtime validation is pending:** run these workflows on the analysis
  machine before interpreting outputs. Synthetic converter regression cases
  can be run there with `python -m unittest discover -s tests -v` in an
  environment containing `cyvcf2`; they do not invoke inference binaries.
- **Relatedness:** `relatedness1` contains IBD segments, so has no enabled
  tools. Other relatedness challenges still need relationship-label/meiosis
  classification from raw estimates; their VCFs are already phased, so
  SHAPEIT5 is not a prerequisite.
- **Reference data:** sweep recombination maps and exon annotations are
  not used by the sweep workflows.
- **Stubs/exclusions:** unsupported tools stay disabled. Reasons and
  implementation requirements live in their configs and module headers.

## Adding a tool

Implement its module, configure installation and arguments in the relevant
challenge configs, then enable it. Existing tool modules are shared across
challenges; new tools also need wiring in the category workflow.

See the [tool catalogue](../docs/existingtools.md) for survey details.
