# GHIST 2026 survey pipeline

Nextflow pipeline for genetic history, sweep detection, and relatedness.
This is a tool-survey scaffold, not a scientifically validated analysis.

## Quick start

Install Nextflow and Conda, then download the challenge data using the
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
| `modules/<category>/<tool>/main.nf` | Shared process implementation per tool |

Tool status is descriptive; `enabled` controls execution:

| Status | Meaning |
|---|---|
| `runnable` | Packaged tool with a documented/standard CLI invocation |
| `runnable-unverified` | Packaged tool; exact invocation still needs verification |
| `runnable-with-caveat` | Runs, but input suitability is limited |
| `excluded` | Not appropriate for the challenge |
| `stub` | Implementation/install work required; disabled by default |

### SweepFinder2 resolution

All four sweep challenges use `args: '-l 1000'`: **1,000 scan positions
across the contig**, not 1,000-base spacing. Change the count in the
challenge config to adjust resolution; this is a coarse survey default.
Bare `-l` is invalid: it omits the required grid count and causes exit 255.
`-lg g` specifies spacing instead.

After changing arguments, use `-resume`. SweepFinder2 reruns conversion,
spectrum estimation, and scanning together because they share one process.

## Known limitations

- **Scientific accuracy is unvalidated.** Check methods against simulations
  with known parameters. Relatedness testing answers are publicly available
  through the [GHIST wiki](https://www.synapse.org/Synapse:syn74803798/wiki/640740);
  final scoring data remains held out.
- **Low-coverage inputs:** BAM archives require Synapse challenge-team
  access. The PSMC/MSMC2 BAM workflows are unverified against real reads.
- **Unverified tools:** Flex-Sweep, fastsimcoal2, and mrpast need their
  exact CLI/data handling checked. See each tool's config and module.
- **Two-population SFS:** fastsimcoal2 currently pools samples rather than
  building a joint 2D spectrum for `secondary_contact` / `split_migration`.
- **Relatedness:** `relatedness1` contains IBD segments, so has no enabled
  tools. Other relatedness challenges still need relationship-label/meiosis
  classification from raw estimates; their VCFs are already phased, so
  SHAPEIT5 is not a prerequisite.
- **Reference data:** recombination maps and sweep exon annotations are
  not wired into modules.
- **Stubs/exclusions:** unsupported tools stay disabled. Reasons and
  implementation requirements live in their configs and module headers.

## Adding a tool

Implement its module, configure installation and arguments in the relevant
challenge configs, then enable it. Existing tool modules are shared across
challenges; new tools also need wiring in the category workflow.

See the [tool catalogue](../docs/existingtools.md) for survey details.
