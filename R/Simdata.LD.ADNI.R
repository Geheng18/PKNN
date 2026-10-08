# ---------------------------------------------------------------------------
# ADNI-panel version of Simdata.LD
#
# Why not UKB: bchQC03 is array data at 7.3 kb/SNP, so a 1000-SNP contiguous
# block spans 6.2 Mb -- a chromosome segment, not a gene.  LD is then so dilute
# that the product kernel K = XX'/m collapses toward constant*I and nothing is
# predictable.  ADNI (the paper's own real data) is dense at 0.07 kb/SNP, so
# 1000 SNPs fit inside ~116 kb of real gene with real LD.
#
# Panel: 24,219 gene-level PLINK filesets, 808 individuals (READ-ONLY on orange).
# Regions are drawn as whole real genes, so "10 genes" means 10 actual genes.
# ---------------------------------------------------------------------------
ADNI.DIR   = "/home/heng.ge/orange/lucienq/adni/QC/adni_gene"
ADNI.INDEX = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/data/adni_gene_index.csv"

#' Draw one replicate from the ADNI gene panel.
#' @param index  data.frame(gene, nSNP) from build_adni_index.R
#' @param n.causal.snp causal SNPs per causal region; NULL = all SNPs in region
Simdata.LD.ADNI = function(index, N, nSNP, nRegion, nRegion.true, simfun,
                           h2 = 0.2, power = 2, n.causal.snp = NULL) {
  suppressMessages(library(BEDMatrix))

  elig = index$gene[index$nSNP >= nSNP]
  if (length(elig) < nRegion)
    stop("only ", length(elig), " genes have >= ", nSNP, " SNPs; need ", nRegion)
  genes = sample(elig, nRegion)

  X.list = vector("list", nRegion); ind = NULL
  for (i in seq_len(nRegion)) {
    bm = BEDMatrix(file.path(ADNI.DIR, paste0(genes[i], ".bed")), simple_names = TRUE)
    if (is.null(ind)) ind = sort(sample(nrow(bm), min(N, nrow(bm))))
    s = sample(ncol(bm) - nSNP + 1, 1)                  # contiguous block
    G = bm[ind, s:(s + nSNP - 1), drop = FALSE]
    if (anyNA(G)) G = apply(G, 2, function(x) { x[is.na(x)] = mean(x, na.rm = TRUE); x })
    keep = apply(G, 2, function(x) length(unique(x)) > 1)
    X.list[[i]] = G[, keep, drop = FALSE]
  }
  N = length(ind)

  y.var = matrix(0, N, N); causal.snp = vector("list", nRegion.true)
  for (i in seq_len(nRegion.true)) {
    Xi = X.list[[i]]
    j  = if (is.null(n.causal.snp)) seq_len(ncol(Xi))
         else sort(sample(ncol(Xi), min(n.causal.snp, ncol(Xi))))
    causal.snp[[i]] = j
    y.var = y.var + findKernel("product", scale(Xi[, j, drop = FALSE], scale = FALSE))
  }
  y.g = MASS::mvrnorm(mu = rep(0, N), Sigma = y.var)
  g   = as.numeric(scale(apply.simfun(y.g, simfun, power)))
  y   = sqrt(h2) * g + sqrt(1 - h2) * rnorm(N)

  info = data.frame(region = seq_len(nRegion), gene = genes,
                    causal = seq_len(nRegion) <= nRegion.true,
                    nSNP.kept = sapply(X.list, ncol),
                    n.causal.snp = c(sapply(causal.snp, length), rep(0L, nRegion - nRegion.true)))
  list(X = X.list, y = y, info = info, causal.snp = causal.snp)
}
