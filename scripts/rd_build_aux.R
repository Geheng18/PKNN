# ---------------------------------------------------------------------------
# Two further inputs the experiments need.
#
# 1. The ELEVEN INDIVIDUAL candidate-gene kernels.  rd_build_kernels.R saves
#    only their sum (K_cand); the interaction basis needs K_1..K_11 separately.
#
# 2. An LD-PRUNED genome-wide genotype matrix for LDpred2, which cannot take
#    20.7M SNPs.  Pruning is done within each gene (MAF >= 0.01, greedy r^2 < 0.2)
#    and the survivors concatenated: between-gene LD is negligible since the
#    genes sit at separate loci, so per-gene pruning is a close approximation to
#    a genome-wide pass and avoids merging 24,219 PLINK filesets.
# ---------------------------------------------------------------------------
suppressMessages({library(BEDMatrix)})
D   = "/home/heng.ge/orange/lucienq/adni/QC/adni_gene"      # READ-ONLY
OUT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/realdata/shared/kernels"
CAND = c("05_141971742","08_038854504","09_134378288","11_005274420","12_053901639",
         "14_094843083","19_045417503","19_045409005","19_045394476","21_036160097","21_048018530")

read.gene = function(g, maf.min=0) {
  X = as.matrix(BEDMatrix(file.path(D, paste0(g,".bed")), simple_names=TRUE))
  if (anyNA(X)) X = apply(X,2,function(x){x[is.na(x)]=mean(x,na.rm=TRUE); x})
  if (maf.min>0) { f=colMeans(X)/2; f=pmin(f,1-f); X=X[, f>=maf.min, drop=FALSE] }
  X[, apply(X,2,function(x) length(unique(x))>1), drop=FALSE]
}

# ---- 1. individual candidate kernels ------------------------------------
K.cand.list = list(); nsnp = integer(0)
for (g in CAND) {
  X = read.gene(g); Xc = scale(X, center=TRUE, scale=FALSE)
  K.cand.list[[g]] = tcrossprod(Xc)/ncol(Xc); nsnp = c(nsnp, ncol(X))
}
cat("individual candidate kernels:", length(K.cand.list), " SNPs:", sum(nsnp), "\n"); flush.console()
save(K.cand.list, nsnp, file=file.path(OUT,"kernels_candidate_individual.Rdata"))

# ---- 2. LD-pruned genome-wide matrix for LDpred2 -------------------------
prune = function(X, r2=0.2) {                    # greedy, keeps first of any correlated pair
  if (ncol(X) < 2) return(seq_len(ncol(X)))
  C = suppressWarnings(cor(X)); C[is.na(C)] = 0
  keep = c(); for (j in seq_len(ncol(X))) if (!length(keep) || all(C[j,keep]^2 < r2)) keep = c(keep,j)
  keep
}
genes = sub("\\.bed$","",list.files(D, pattern="\\.bed$"))
cols = list(); labs = character(0); t0 = Sys.time()
for (i in seq_along(genes)) {
  X = read.gene(genes[i], maf.min=0.01)
  if (!ncol(X)) next
  k = prune(X)
  cols[[length(cols)+1]] = X[, k, drop=FALSE]
  labs = c(labs, paste0(genes[i],":",k))
  if (i %% 2000 == 0) { cat(sprintf("  %5d/%d genes, %d SNPs kept (%.1f min)\n", i, length(genes),
      sum(sapply(cols,ncol)), as.numeric(difftime(Sys.time(),t0,units="mins")))); flush.console() }
}
G.pruned = do.call(cbind, cols); colnames(G.pruned) = labs
gene.of  = sub(":.*","",labs)
cat("pruned genome-wide matrix:", dim(G.pruned), "\n")
save(G.pruned, gene.of, file=file.path(OUT,"pruned_genomewide.Rdata"))
cat("DONE in", round(as.numeric(difftime(Sys.time(),t0,units="mins")),1), "min\n")
