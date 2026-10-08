# ---------------------------------------------------------------------------
# Interaction-effect generator, identical to Simdata.LD.R except that the
# genetic value is driven by explicit pairwise SNP products instead of a
# nonlinear link function.  Same UKB panel, same contiguous one-region-per-
# chromosome LD blocks, same explicit heritability targeting, same output layout.
#
# Ported from within.region()/outside.region() in KNN_Data/R/Simdata.R:
#   within.region   pairs of SNPs drawn inside each causal gene
#   outside.region  pairs drawn across the pooled causal genes
#   threshold=TRUE  negative products clipped to zero (a threshold effect);
#   threshold=FALSE the raw product is kept (a multiplicative effect)
#
# NOTE ON NAMING.  In the surviving Simdata.R the dispatch is crossed --
#   switch(inter, within.region = {outside.region(...)},
#                 outside.region = {within.region(...)})
# so 2023 directories may hold the opposite construction to their name.  That
# code path was not the one used for the published data (Simdata.interaction,
# now lost), so the crossing cannot be confirmed.  Here the name matches what is
# actually built: within.region.* really is within-gene pairing.
# ---------------------------------------------------------------------------

# pairwise products inside each causal gene
inter.within = function(X.list, nRegion.true, threshold) {
  lapply(seq_len(nRegion.true), function(i) {
    X = scale(X.list[[i]], scale = FALSE)
    m = ncol(X); np = m %/% 2
    if (np < 1) stop("gene too small to pair")
    idx = matrix(sample(m, 2 * np), np, 2)
    Z = X[, idx[,1], drop = FALSE] * X[, idx[,2], drop = FALSE]
    if (threshold) Z[Z < 0] = 0
    Z })
}
# pairwise products across the pooled causal genes
inter.outside = function(X.list, nRegion.true, threshold) {
  P = do.call(cbind, lapply(X.list[seq_len(nRegion.true)], function(G) scale(G, scale = FALSE)))
  m = ncol(P); np = m %/% 2
  idx = matrix(sample(m, 2 * np), np, 2)
  Z = P[, idx[,1], drop = FALSE] * P[, idx[,2], drop = FALSE]
  if (threshold) Z[Z < 0] = 0
  per = np %/% nRegion.true
  lapply(seq_len(nRegion.true), function(i) Z[, ((i-1)*per + 1):(i*per), drop = FALSE])
}

#' @param scenario one of within.region.threshold.TRUE / .FALSE /
#'                 outside.region.threshold.TRUE / .FALSE
#' @param n.causal.snp number of causal SNPs inside each causal region that may
#'   enter an interaction pair.  NULL = every SNP in the region is eligible (the
#'   original design).  A fixed count gives the WITHIN-REGION sparsity Reviewer 2
#'   (comment 1) asked for: the model still receives all nSNP columns, but only
#'   this many of them drive the phenotype, so the rest enter as noise.
Simdata.LD.INT = function(bm, bim, N, nSNP, nRegion, nRegion.true, scenario, h2 = 0.2,
                          n.causal.snp = NULL) {
  parts     = strsplit(scenario, "\\.threshold\\.")[[1]]
  inter     = parts[1]; threshold = as.logical(parts[2])
  stopifnot(inter %in% c("within.region","outside.region"), !is.na(threshold))

  ind  = sort(sample(nrow(bm), N))
  cols = pick.regions(bim, nSNP, nRegion)
  X.list = vector("list", nRegion)
  for (i in seq_len(nRegion)) {
    G = bm[ind, cols[[i]], drop = FALSE]
    keep = apply(G, 2, function(x) length(unique(x[!is.na(x)])) > 1)
    if (any(!keep)) G = G[, keep, drop = FALSE]
    if (anyNA(G)) G = apply(G, 2, function(x) { x[is.na(x)] = mean(x, na.rm = TRUE); x })
    X.list[[i]] = G
  }

  # Within-region sparsity: restrict which SNPs are eligible to be paired.  The
  # subset drives the phenotype; X.list (all nSNP columns) is what the model
  # sees, so the unselected SNPs act as noise inside a causal region.
  causal.snp = vector("list", nRegion.true)
  C.list = X.list
  if (!is.null(n.causal.snp)) {
    for (i in seq_len(nRegion.true)) {
      j = sort(sample(ncol(X.list[[i]]), min(n.causal.snp, ncol(X.list[[i]]))))
      causal.snp[[i]] = j
      C.list[[i]] = X.list[[i]][, j, drop = FALSE]
    }
  } else {
    for (i in seq_len(nRegion.true)) causal.snp[[i]] = seq_len(ncol(X.list[[i]]))
  }

  Z.list = if (inter == "within.region") inter.within(C.list, nRegion.true, threshold)
           else                          inter.outside(C.list, nRegion.true, threshold)
  y.var = matrix(0, N, N)
  for (Z in Z.list) y.var = y.var + findKernel("product", scale(Z, scale = FALSE))

  y.g = MASS::mvrnorm(mu = rep(0, N), Sigma = y.var)
  g   = as.numeric(scale(y.g))
  if (!is.finite(sd(g)) || sd(g) == 0) stop("degenerate signal for ", scenario)
  y = sqrt(h2) * g + sqrt(1 - h2) * rnorm(N)

  info = data.frame(region = seq_len(nRegion),
                    causal = seq_len(nRegion) <= nRegion.true,
                    chr = sapply(cols, function(z) bim$V1[z[1]]),
                    start.bp = sapply(cols, function(z) bim$V4[z[1]]),
                    end.bp = sapply(cols, function(z) bim$V4[z[length(z)]]),
                    nSNP.kept = sapply(X.list, ncol),
                    n.inter.pairs = c(sapply(Z.list, ncol), rep(0L, nRegion - nRegion.true)),
                    n.causal.snp = c(sapply(causal.snp, length), rep(0L, nRegion - nRegion.true)))
  list(X = X.list, y = y, info = info)
}
