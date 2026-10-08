# PKNN — A Penalized Kernel Neural Network for Complex Genetic Risk Prediction

Implementation, simulation studies and results for the paper.

## What a kernel neural network is here

Genotypes are grouped into biologically defined regions, each contributing one
input kernel. A hidden layer combines those kernels through random weights, so
the induced covariance of the network output is a linear combination of kernel
products and the network is equivalent to a linear mixed model whose variance
components are the combination weights. Those components are estimated and
selected by penalized MINQUE. With ten regions the interaction basis has 67
components: an identity term, ten main-effect kernels, and the Hadamard products
between pairs of regions.

Under a polynomial kernel and an L2 penalty the MINQUE solution is closed form,
which is where the computational advantage over iterative REML comes from. L1 and
elastic-net penalties are supported through numerical optimisation.

## Layout

```
R/            method implementation
scripts/      simulation generators, analysis drivers, figure scripts, SLURM jobs
results/      per-replicate evaluation output for every figure in the paper
figures/      the figures themselves, plus figures.tex with captions
data/         synthetic genotype generator; see data/README.md for why no real data
```

### Method (`R/`)

| File | Purpose |
|---|---|
| `KNN2LMM.R` | expands a kernel neural network into its LMM variance components |
| `MINQUE.selection.R` | closed-form MINQUE under an L2 penalty |
| `MINQUE.selection.numerical.R` | numerical MINQUE under L1 / elastic-net penalties |
| `MINQUE.selection.cv.R` | cross-validated choice of the penalty parameter |
| `MINQUE.predict.R` | prediction for held-out individuals |
| `GMMLasso.R`, `ChooseLambda.R` | the GMM comparator |
| `KernelPool.R`, `LinearKernel.R`, `ScaleKernel.R` | kernel construction |
| `Simdata.LD.R` | simulation under nonlinear architectures, preserving real LD |
| `Simdata.LD.INT.R` | simulation under pairwise interaction architectures |

### Simulation studies (`scripts/`)

| Study | Generator | Analysis |
|---|---|---|
| Main grid: 9 architectures x 4 SNP densities x 5 causal rates, two heritabilities | `gen_sim_data.R`, `gen_sim_int.R` | `analyze_sim.R`, `analyze_extra.R` |
| Number of candidate genes: 2 to 20 regions supplied, 2 causal | `gen_sim_ngene.R`, `gen_sim_ngene_int.R` | `analyze_ngene.R`, `analyze_ngene_extra.R` |
| Within-region sparsity: 5, 10, 15 or 20 causal SNPs per causal region | `gen_sim_sparse.R` | `analyze_sparse.R`, `analyze_sparse_extra.R` |
| Sample size: N = 500 to 5,000 (10,000 for the linear case) | `sim_nvary.R` | same file |
| Runtime and memory: N = 500 to 10,000 x 20 to 10,000 SNPs per gene | `benchmark_NP.R` | same file |

Figures are produced by `make_figures.R`, `make_sim_panels.R`,
`make_ngene_figure.R`, `make_sparse_figure.R` and `make_nvary_figure.R`.
`make_realdata_tables.R` produces the application tables.

### Comparators

Closed-form MINQUE ridge, numerical MINQUE lasso and elastic net, GMM, gBLUP,
LDpred2 (`bigsnpr`), hierNet, and a covariates-only baseline. glinternet was
also evaluated; it is reported in the paper's text rather than its figures
because it returned a usable fit in only 74% of replicates even at the smallest
region size. SHIM was requested in review but has no maintained public
implementation and could not be run.

### Application (`scripts/rd_*.R`)

The ADNI pipeline is included as code. ADNI data is not redistributable, so these
scripts will not run without your own access and will need their input paths
changed. `rd_setup.R` builds the phenotypes and splits, `rd_build_kernels.R` the
candidate-gene and rest-of-genome kernels, `rd_analyze.R` runs the comparison,
and `rd_lambda_cv.R` tunes the ridge parameter inside each training set.

## Reproducing

Scripts are written for SLURM; the `.sbatch` files record the resources each
stage needed. To run a single simulation cell interactively, set
`SLURM_ARRAY_TASK_ID` and call the script directly:

```sh
export SLURM_ARRAY_TASK_ID=1
Rscript scripts/gen_sim_data.R      # then
Rscript scripts/analyze_sim.R
```

Paths at the top of each script point at the authors' cluster and need editing.
R 4.1 was used, with `MASS`, `glmnet`, `matrixcalc`, `Matrix`, `rrBLUP`,
`bigsnpr`, `bigsparser`, `hierNet`, `BEDMatrix`, `data.table`, `ggplot2` and
`ggpubr`.

## Results

`results/` holds per-replicate output — one row per replicate, method and cell,
with mean squared error, correlation and wall-clock time — so every figure can be
re-derived without recomputation. The aggregated summaries used directly by the
figures are `ngene_summary.csv`, `sparse_summary.csv`, `nvary_summary.csv`,
`figure_validity.csv` and `panel_validity_P20.csv`.

A caution when reading the results: GMM frequently returns a degenerate fit
whose genetic component has been shrunk to a constant. Its means are therefore
taken over the subset of replicates that did not collapse — a third or fewer in
several settings — and the validity files record the usable count for every cell.

## Data

No genotype data is distributed. See `data/README.md`.
