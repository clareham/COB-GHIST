# GHIST 2026 survey pipeline

An early-survey Nextflow pipeline across the [GHIST 2026](../docs/existingtools.md)
challenge set, using the data dump in `data-GHIST2026/` (a sibling folder
to this repo). The goal of this scaffold is breadth: get a real,
installable tool wired up against every challenge it plausibly applies
to, and be explicit about the rest (data/tool mismatches, source-only
builds, unverified CLI usage) rather than quietly skipping them.

**Only `*.testing.*` / `*_testing.*` files are ever read.** `*.final.*`
files are the held-out GHIST scoring sets and must never be consumed by
this pipeline - enforced centrally by `GhistUtils.testingFile()` (see
below), not by convention alone.

## Config hierarchy

The whole pipeline is organized as three nested levels, mirrored 1:1 in
`conf/`:

```
conf/
├── base.config                        # shared process defaults
├── genetic_history/
│   ├── genetic_history.config         # level 1: category
│   ├── bottleneck.config              # level 2: challenge (+ level 3 `tools {}` block)
│   ├── bottleneck_low_coverage.config
│   ├── growth.config
│   ├── secondary_contact.config
│   ├── split_migration.config
│   └── split_isolation_low_coverage.config
├── sweep_detection/
│   ├── sweep_detection.config         # level 1: category
│   ├── onesweep.config                # level 2 + 3
│   ├── onesweep_bgs.config
│   ├── softsweep.config
│   └── multisweep_bgs.config
└── relatedness/
    ├── relatedness.config             # level 1: category
    ├── relatedness1.config            # level 2 + 3 (empty tools - see below)
    ├── relatedness2.config
    ├── relatedness3.config
    └── relatedness4.config
```

- **Level 1 (category)**: `genetic_history` / `sweep_detection` /
  `relatedness`. A single `params.run_<category>` boolean (set in
  `nextflow.config`) turns the whole category on/off.
- **Level 2 (challenge)**: one config per GHIST scenario, each with its
  own `enabled` flag, input file name(s), and any challenge-specific
  reference data (accessible-region BED, coding-region BED, submission
  schema path).
- **Level 3 (tool)**: a `tools { <tool_name> { enabled, status, conda,
  args, ... } }` block nested inside each challenge config - one entry
  per tool workflow from `docs/existingtools.md` that's a plausible fit
  for that challenge's input format.

Code is **not** duplicated per challenge: each tool has exactly one
process, under `modules/<category>/<tool>/main.nf`, reused across every
challenge in that category. The three `workflows/<category>.nf` files do
the fan-out - for every enabled (challenge, tool) pair in the config,
they route that challenge's testing file into that tool's process. A
tool with zero enabled challenges in a given run just gets an empty
input channel and never executes (no extra conditional logic needed) -
this is how every `stub`/`excluded` entry stays inert by default.

```
conf/<category>/<challenge>.config   →  params.<category>.challenges.<challenge>.tools.<tool>
                                                        │
workflows/<category>.nf  (fans out, builds (meta, file) per enabled pair)
                                                        │
modules/<category>/<tool>/main.nf    (one process per tool, shared across challenges)
```

### Why the testing/final guard lives in `lib/`, not `conf/`

Nextflow's config files use a restricted declarative DSL (no custom
`def` functions/classes) - only `lib/*.groovy` classes and the
`workflows/*.nf` files can run arbitrary Groovy. So challenge configs
only declare a bare filename (e.g. `testing_vcf_name: 'bottleneck.testing.vcf'`),
and each `workflows/<category>.nf` resolves it via
`GhistUtils.testingFile(params.data_dir, name)` (`lib/GhistUtils.groovy`),
which throws if the name doesn't look like a testing file. This was
confirmed by a real parse/run error during development, not a guess - see
`lib/GhistUtils.groovy`'s docstring for the exact rule (it also handles
`relatedness2_testing.vcf`'s underscore instead of the usual dot).

## Tool status legend

Every `tools {}` entry carries a `status`, set deliberately rather than
uniformly "enabled: true":

| status | meaning |
|---|---|
| `runnable` | real conda/pip package, CLI usage grounded in documented/standard conventions |
| `runnable-unverified` | real package, but the exact CLI invocation for this specific use case (data format, binary name, flags) is a best-effort guess, not cross-checked against the tool's own docs/manual |
| `runnable-with-caveat` | runs mechanically, but the *input data itself* doesn't give the tool a fair shot (documented per-entry) |
| `excluded` | package may exist, but this challenge's data doesn't suit the tool's intended use case |
| `stub` | no conda/PyPI package found; needs a manual from-source install before it can run at all |

None of this has been validated for scientific correctness end-to-end
(that needs a round-trip against data with known simulation parameters -
see Limitations). What *has* been verified by actually running this
pipeline in this environment:
- The full config hierarchy parses cleanly (`nextflow config .`).
- Channel fan-out produces exactly the expected task count per tool
  across all three categories (e.g. `PSMC`: 6 tasks, one per applicable
  genetic_history challenge; `SWEEPFINDER2`/`FLEXSWEEP`: 4 tasks each
  across the four sweep challenges; `KING_PLINK`/`NGSRELATE`/
  `ADMIXTURE`/`SHAPEIT5`: 3 tasks each, correctly skipping
  `relatedness1`; every `stub`/`excluded` tool: 0 tasks).
- `PSMC` on `bottleneck_low_coverage` runs for real end to end via
  `-profile conda` against the actual testing FASTA and produces a real
  `.psmc` file - which (as documented) comes back with `sum_n:0`
  (zero heterozygous windows detected), i.e. exactly the flat/uninformative
  curve the caveat predicts, since that challenge only hands out a single
  haploid sequence with no het signal.

## Running it

```console
# Parse/validate the full config hierarchy without running anything
nextflow config .

# Smoke-test a single category with conda-managed per-tool environments
nextflow run . -profile conda --run_sweep_detection false --run_relatedness false

# Full survey run
nextflow run . -profile conda
```

Toggle scope at any level without touching code:
- Whole category: `--run_genetic_history false`
- One challenge: edit `enabled` in its `conf/<category>/<challenge>.config`
- One tool for one challenge: edit `enabled` in that challenge's `tools {}` block

## Known limitations (read before trusting any output)

- **No ground truth is available locally** (by design - `*.final.*` is
  held out), so there's no way to score inference accuracy here. Validate
  methods first against your own `msprime`/`stdpopsim` simulations with
  known parameters, then apply to the real testing data.
- **No popmap for two-population challenges** (`secondary_contact`,
  `split_migration`): `fastsimcoal2`'s SFS is pooled into a single
  population rather than a true joint 2D SFS (see `bin/vcf_to_sfs.py`).
- **`*_low_coverage` challenges ship a single haploid consensus FASTA**,
  not reads/BAM or a diploid pair - `psmc` runs on it but produces a
  structural placeholder, not a validated inference (verified above);
  `msmc2` is excluded outright for lacking any het signal to use.
- **`relatedness1`'s testing data is already pairwise IBD segments**, not
  genotypes - none of the catalogued relatedness tools take that as
  primary input, so its `tools {}` map is intentionally empty pending a
  clearer idea of the intended task/tool.
- **RAiSD, SweeD, OmegaPlus, Stairway Plot 2, fastGLOBETROTTER, ArchIE2,
  READ/KIN, IcMLkin, GRUPS-rs, PHLASH** are `stub`s: no conda/PyPI
  package was found, and (per explicit direction) compiling tools like
  RAiSD from source wasn't worth the time for this pass. Each stub
  module documents exactly what a real implementation needs.
- `fastsimcoal2` and `mrpast`'s exact CLI/model-file syntax is a
  best-effort guess (`status: 'runnable-unverified'`) - the bioconda
  `fastsimcoal2` binary is version-tagged (resolved dynamically in the
  module) and `mrpast` is ARG-based, so its VCF input is first run
  through `tsinfer` inline.

## Extending

To turn a `stub` into something real: install the tool (conda/pip/source
build), fill in its module's `script:` block, set `conda`/`container` in
the relevant challenge config(s), and flip `enabled: true`. No workflow
or channel-wiring changes are needed - the fan-out in
`workflows/<category>.nf` already routes any enabled (challenge, tool)
pair automatically.
