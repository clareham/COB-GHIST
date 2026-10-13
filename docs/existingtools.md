# Notes on Existing Tools

## Genetic History Challenges

Common methods include: `Ancestral Recombination Graphs (ARG)`, `Hidden Markov Model Coalescence`, `Site Frequency Spectrum`, `Linkage Disequilibrium`.

### PHLASH

- Python API
- Population History Learning by Averaging Sampled Histories

### mrpast

- composite tool utilizing ARGs, HMMs, and History Simulation testing
- Can integrate with dedicated ARG tools or infer parameters

### MSMC2 / PSMC

- Pairwise / Multiple Sequentially Markovian Coalescent
- Unclear what expected inputs are, noted that its designed to work with WGS

### Stairway Plot 2

- SFS, CLI-only

### fastsimcoal2

- SFS + Markov

### fastGLOBETROTTER

- linkage-disequilibrium

### ArchIE2

- Identification of introgressed DNA (i.e. from another species)
- Came up in search, may be tangentially useful?

## Sweep Detection Challenges

### SweepFinder / SweepFinder2

- original / modernized versions (2016)
- Utilizes SFS
- Includes controls for variable mutation rate and background noise
- One of earliest sweep detection tools, now faster in modern implementation

### SweeD

- 2005
- Derived from original SweepFinder
- Much faster, fewer controls

### OmegaPlus

- linkage-disequilibrium
- computes composite omega statistic across sliding window
- site-specific

### RAiSD

- composite score from VMR, SFS, and LD Omega
- maintenance is not great, appears need to compile from source
- Depends on GNU scientific library native install

### Flex-Sweep 2

- CNN with 1-d kernel
- stacks 11 summary stats into matrix across sliding window sampling
- Domain-Adaptive Neural Network approach allows fine-tuning before use

### Timesweeper

- Similar approach to Flex-Sweep 2, designed for time-series data.

## Relatedness Estimation Challenges

### KING & PLINK

- Kinship-based Inference for GWAS
- Estimates pairwise kinship coefficients

### READ & KIN

- Modernized equivalents utilizing HMM coalescence.

### NgsRelate

- Optimized for low-coverage ngs data
- Returns likelihoods instead of binary calls

### IcMLkin

- low-coverage maximum-likelihood kinship estimation
- requires raw reads, attempts to correct for low-coverage and ancient DNA with poorly related reference genomes
- Linkage-Disequilibrium

### GRUPS-rs

- pedigree simulations for background estimation
- rust

### SHAPEIT5

- Available in Python & R

### ADMIXTURE

- Modernized implementation of STRUCTURE algo
- considered current SOTA
- Operates with or without priors (i.e. expected lineage information)
- Utilizes unsupervised clustering if not provided priors
