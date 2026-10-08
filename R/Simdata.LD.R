# ---------------------------------------------------------------------------
# Simdata.LD -- simulation-data generator preserving TRUE LD
#
# Rewrite of the lost `Simdata.nonlinear` used by data_generate1.R (blocks at
# lines 104-142), following the surviving `Simdata` / data_generate_test_10_region.R
# template but with two changes required by the reviewers:
#
#   1. TRUE LD.  The original drew SNPs at random from within a ~500-SNP window
#      (calibrated from the saved data: within-gene mean r^2 = 0.0129 vs 0.0044
#      for genome-wide random and 0.0446 for fully contiguous).  Here each
#      region is a *contiguous* run of nSNP SNPs, and the 10 regions are placed
#      on 10 distinct chromosomes so between-region LD is nil.
#
#   2. EXPLICIT HERITABILITY.  The original mixed signal and noise through
#      raneff.scale / noise.scale, which had to be retuned per link function
#      because cosh/hyperbola/rickercurve change the variance.  Here the
#      transformed signal is standardised and mixed at a target h2, so effect
#      size is comparable across link functions and directly reportable.
#
# Genotypes: UK Biobank bchQC03 (61,463 individuals x 341,545 SNPs), read
# lazily via BEDMatrix -- no plink subsetting calls needed.
# ---------------------------------------------------------------------------

# --- nonlinear link functions (identical to KNN_Data/R/Simdata.R) -----------
r.softplus  = function(x) log(1 + exp(x))
f.power     = function(x, power = 2) x^power
f.cosh      = function(x) cosh(x)
f.hyperbola = function(x, power = 2) 20 * (r.softplus(x^power) / (1 + r.softplus(x^power)))
f.rickercurve = function(x, power = 2) 60 * (r.softplus(x^power) * exp(-r.softplus(x^power)))

apply.simfun = function(y, simfun, power = 2) {
  switch(simfun,
    linear      = y,
    square      = f.power(y, power),
    power       = f.power(y, power),
    cosh        = f.cosh(y),
    hyperbola   = f.hyperbola(y, power),
    rickercurve = f.rickercurve(y, power),
    stop("unknown simfun: ", simfun))
}

#' Pick nRegion contiguous SNP blocks, one per chromosome
#'
#' @param bim   data.frame with columns V1 (chr) and V4 (bp)
#' @param nSNP  SNPs per region
#' @param nRegion number of regions
#' @return list of integer vectors of column indices into the .bed
pick.regions = function(bim, nSNP, nRegion) {
  chrs = unique(bim$V1)
  # chromosomes with enough SNPs to hold a block
  ok = chrs[sapply(chrs, function(c) sum(bim$V1 == c) >= nSNP + 1)]
  if (length(ok) < nRegion)
    stop("only ", length(ok), " chromosomes can hold ", nSNP, " SNPs; need ", nRegion)
  use = sample(ok, nRegion)
  lapply(use, function(c) {
    idx = which(bim$V1 == c)
    s   = sample(length(idx) - nSNP + 1, 1)
    idx[s:(s + nSNP - 1)]
  })
}

#' Generate one replicate
#'
#' @param bm       BEDMatrix handle
#' @param bim      bim table
#' @param N        individuals
#' @param nSNP     SNPs per region
#' @param nRegion  total regions
#' @param nRegion.true number of causal regions (causal rate = nRegion.true/nRegion)
#' @param simfun   link function name
#' @param h2       target heritability on the transformed scale
#' @return list(X = list of genotype matrices, y = phenotype, info = region metadata)
#' @param n.causal.snp number of causal SNPs inside each causal region.
#'   NULL = every SNP in the region contributes (the original design).
#'   A fixed count makes larger nSNP add *noise* SNPs rather than diluting the
#'   signal.  This is required because the product kernel K = XX'/m has
#'   off-diagonal SD proportional to 1/sqrt(m): as m grows K -> constant * I, the
#'   simulated genetic value becomes i.i.d. noise, and nothing is predictable.
#'   It also gives the within-region sparsity Reviewer 2 (comment 1) asked for.
Simdata.LD = function(bm, bim, N, nSNP, nRegion, nRegion.true, simfun,
                      h2 = 0.2, power = 2, n.causal.snp = NULL) {

  ind  = sort(sample(nrow(bm), N))
  cols = pick.regions(bim, nSNP, nRegion)

  X.list = vector("list", nRegion)
  for (i in seq_len(nRegion)) {
    G = bm[ind, cols[[i]], drop = FALSE]
    # monomorphic SNPs carry no signal and break the kernel scaling
    keep = apply(G, 2, function(x) length(unique(x[!is.na(x)])) > 1)
    if (any(!keep)) G = G[, keep, drop = FALSE]
    if (anyNA(G)) G = apply(G, 2, function(x) { x[is.na(x)] = mean(x, na.rm = TRUE); x })
    X.list[[i]] = G
  }

  # genetic value: MVN with covariance = sum of causal-region product kernels
  y.var = matrix(0, N, N)
  causal.snp = vector("list", nRegion.true)
  for (i in seq_len(nRegion.true)) {
    Xi = X.list[[i]]
    if (is.null(n.causal.snp)) {
      j = seq_len(ncol(Xi))
    } else {
      j = sort(sample(ncol(Xi), min(n.causal.snp, ncol(Xi))))
    }
    causal.snp[[i]] = j
    Xc = scale(Xi[, j, drop = FALSE], scale = FALSE)
    y.var = y.var + findKernel("product", Xc)
  }
  y.g = MASS::mvrnorm(mu = rep(0, N), Sigma = y.var)

  # nonlinear transform, then mix signal and noise at the target heritability
  g = apply.simfun(y.g, simfun, power)
  g = as.numeric(scale(g))                       # unit variance
  if (!is.finite(sd(g)) || sd(g) == 0) stop("degenerate signal for simfun=", simfun)
  y = sqrt(h2) * g + sqrt(1 - h2) * rnorm(N)

  info = data.frame(
    region = seq_len(nRegion),
    causal = seq_len(nRegion) <= nRegion.true,
    chr    = sapply(cols, function(z) bim$V1[z[1]]),
    start.bp = sapply(cols, function(z) bim$V4[z[1]]),
    end.bp   = sapply(cols, function(z) bim$V4[z[length(z)]]),
    nSNP.kept = sapply(X.list, ncol))
  info$n.causal.snp = c(sapply(causal.snp, length), rep(0L, nRegion - nRegion.true))

  list(X = X.list, y = y, info = info, causal.snp = causal.snp)
}
