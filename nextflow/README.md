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
  across all three categories (e.g. `PSMC`: 4 VCF-based tasks, correctly
  skipping the two BAM-based `*_low_coverage` challenges when their
  archives aren't present locally; `SWEEPFINDER2`/`FLEXSWEEP`: 4 tasks
  each across the four sweep challenges; `KING_PLINK`/`NGSRELATE`/
  `ADMIXTURE`/`SHAPEIT5`: 3 tasks each, correctly skipping
  `relatedness1`; every `stub`/`excluded` tool: 0 tasks).
- `PSMC` on a real `-profile conda` run produced a genuine `.psmc` file
  end to end - that first real run (before the wiki cross-check below)
  used `bottleneck_low_coverage`'s bare reference FASTA directly and came
  back with `sum_n:0` (zero heterozygous windows), which is exactly what
  you'd expect from a single haploid sequence with no het signal to find.
  That result is *why* the data format was cross-checked against the
  official GHIST wiki (see next section) - it turned out the real
  challenge data is low-coverage BAMs, not a bare FASTA, and the module
  has since been rewritten accordingly (untested against real reads, see
  Limitations).
- That rewrite surfaced a real Nextflow bug during testing: passing the
  same file object to two separate `path` inputs caused a work-directory
  staging name collision, fixed via `stageAs` (see
  `modules/genetic_history/psmc/main.nf`).

## Data validation against the official GHIST wiki

The [GHIST 2026 Synapse wiki](https://www.synapse.org/Synapse:syn74803798/wiki/640740)
documents the exact data format/metadata for every challenge. Cross-
referencing every file ID it links against `manifest.csv` found our local
data drop has **41 of 55** referenced files; the other 14 fall into three
buckets:

- **High impact, since fixed in config/code but not yet runnable**: the
  `*_low_coverage` challenges' real data is a `tar.gz` archive of 20
  diploid individuals' aligned BAMs (~2-4x coverage) + a `popfile.txt`,
  not the bare consensus FASTA this pipeline originally assumed (see
  below). The archives (`syn77432049`, `syn77432048`, `syn77432188`,
  `syn77432194`) require Synapse challenge-team authentication this
  environment doesn't have - `psmc`/`msmc2` have been rewritten for the
  real BAM format and will run automatically once you download the
  archives into `data_dir` yourself, but that rewrite is unverified
  against real reads in the meantime.
- **Medium impact, genuinely missing, pipeline still runs without them**:
  recombination maps for `growth`/`split_migration` (`syn75318370`/
  `syn75318369`) and for sweep detection chr7/chr9 (`syn77026803`/
  `syn77026807`), plus sweep detection's exon/Havana annotations
  (`syn77027375`/`syn77027379`). Not wired into any module yet.
- **Not actually a gap - an opportunity**: unlike every other category,
  the relatedness challenges' *testing* ground truth is itself a public,
  downloadable file (`syn75202646`, `syn75158967`, `syn75203704`,
  `syn75203706`), not just inline-revealed numbers. Worth downloading to
  validate `king_plink`/`ngsrelate`/`admixture` output - something no
  other category allows.

Two correctness bugs were also caught this way and fixed:
1. `growth.config` was missing `coding_regions_bed` - the wiki confirms
   "Growth" is explicitly "Growth **with Background Selection**", using
   the same BED as `split_migration`.
2. (Not yet fixed - flagged for a future pass) relatedness VCFs are
   documented as "phased and error-free", so `shapeit5`'s "phasing
   prerequisite" framing is inaccurate for this data; and the actual
   expected *output* for relatedness2/3/4 is a relationship label/meiosis
   count per pair (siblings, avuncular, Nth cousin, ...), which
   `king_plink`/`ngsrelate`'s raw kinship coefficients don't directly
   provide - there's no classification step from coefficient -> label
   yet.

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
  known parameters, then apply to the real testing data. (Exception:
  the relatedness challenges' testing answers are public - see above.)
- **No popmap for two-population challenges** (`secondary_contact`,
  `split_migration`): `fastsimcoal2`'s SFS is pooled into a single
  population rather than a true joint 2D SFS (see `bin/vcf_to_sfs.py`).
- **The `*_low_coverage` challenges' BAM archives aren't downloaded**
  (require Synapse challenge-team auth - see above). `psmc`/`msmc2` are
  wired up for the real BAM+popfile format (`bin/bam_sample_to_psmcfa.sh`,
  `bin/bam_archive_to_merged_vcf.sh`, `bin/bam_to_depth_mask.py`) and will
  run automatically once the archives are in `data_dir`, but that whole
  path is `status: 'runnable-unverified'` - none of it has executed
  against real reads yet.
- **`relatedness1`'s testing data is already pairwise IBD segments**, not
  genotypes - none of the catalogued relatedness tools take that as
  primary input, so its `tools {}` map is intentionally empty pending a
  clearer idea of the intended task/tool.
- **Relatedness output gap**: the challenges want a relationship
  label/meiosis count per pair, not just a raw kinship coefficient -
  `king_plink`/`ngsrelate`/`admixture` only provide the latter; no
  classification step exists yet.
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
