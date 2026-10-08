# Data

## What is here

Nothing but a synthetic genotype generator. No genotype data is distributed with
this repository.

## Why

The simulations draw genotypes from **UK Biobank** and the application uses
**ADNI**. Neither may be redistributed, so neither is included.

This restriction extends further than it may first appear. The simulation
"data" produced by `scripts/gen_sim_*.R` is only partly synthetic: the
**phenotypes** are simulated, but each `Sim_<rep>_<region>_<n>.Rdata` file holds
a block of **real UK Biobank genotype dosages** for the sampled individuals.
Those files are therefore UK Biobank data and are not shared. The phenotypes
alone would be of no use without them.

## What to do instead

**With UK Biobank access.** The simulations are exactly reproducible. Every
generator seeds deterministically from its array index and replicate number
(see `set.seed` in each `gen_sim_*.R`), so pointing `UKB` at the same genotype
file regenerates the published datasets bit for bit. The file used was a
quality-controlled subset of 61,463 individuals and 341,545 SNPs.

**Without UK Biobank access.** `synthetic_genotypes.R` builds a panel with the
structure the pipeline relies on: contiguous SNP windows carrying linkage
disequilibrium, one window per chromosome. Replace the two lines that open the
real data in any `gen_sim_*.R`:

```r
source("data/synthetic_genotypes.R")
syn = synthetic_panel(nInd = 5000, nChr = 20, nSNPperChr = 2000)
bm = syn$bm; bim = syn$bim
```

The pipeline then runs end to end. Results will not match the published numbers,
because the genotypes are not the same data; the default `rho = 0.85` gives a
mean within-window pairwise r^2 near the 0.045 measured in the UK Biobank
windows used in the paper, against about 0.004 for SNPs drawn at random.
Check with `window_r2(syn$bm, syn$bim)`.

## Results are included

Every number behind every figure is in `results/` as per-replicate evaluation
output, so the analysis can be verified, re-aggregated or re-plotted without
re-running anything.
