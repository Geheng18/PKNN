# ---------------------------------------------------------------------------
# Accumulate genome-scale kernels for the ADNI real-data experiments.
#
# The rest of the genome is 20,722,732 SNPs across 24,208 genes.  As a matrix
# that is 808 x 20.7M = 125 GB, but the product kernel is additive over SNPs:
#       K = sum_g X_g X_g' / M
# so it is accumulated gene by gene into a single 808 x 808 matrix (5 MB).
# Every gene fileset carries the same 808 individuals in the same order
# (verified: 400/400 .fam files byte-identical), so no ID matching is needed.
#
# Genotypes are centred per SNP across all 808 individuals -- one global kernel
# reused by every split.  Recomputing per split would mean 600 passes over
# 20.7M SNPs.  This matches standard GRM practice.
#
# Produces, all 808 x 808:
#   K_cand      11 candidate genes          (10,011 SNPs)
#   K_rest      24,208 remaining genes      (20.7M SNPs)
#   K_rest_maf01  same, MAF >= 0.01
#   K_all       all 24,219 genes            (whole-genome GRM for BLUP)
# ---------------------------------------------------------------------------
suppressMessages({library(BEDMatrix)})
D  = "/home/heng.ge/orange/lucienq/adni/QC/adni_gene"      # READ-ONLY
OUT= "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/realdata/shared/kernels"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

CAND = c("05_141971742","08_038854504","09_134378288","11_005274420","12_053901639",
         "14_094843083","19_045417503","19_045409005","19_045394476","21_036160097","21_048018530")
genes = sub("\\.bed$","",list.files(D, pattern="\\.bed$"))
rest  = setdiff(genes, CAND)
cat("genes:", length(genes), " candidate:", length(CAND), " rest:", length(rest), "\n"); flush.console()

fam = read.table(file.path(D, paste0(CAND[1],".fam")), stringsAsFactors=FALSE)
n   = nrow(fam); ids = fam$V2
acc = function() list(K=matrix(0,n,n), m=0L)

add = function(a, g, maf.min=0) {
  bm = BEDMatrix(file.path(D, paste0(g,".bed")), simple_names=TRUE)
  X  = as.matrix(bm)
  if (anyNA(X)) X = apply(X, 2, function(x){ x[is.na(x)] = mean(x, na.rm=TRUE); x })
  if (maf.min > 0) {
    f = colMeans(X)/2; f = pmin(f, 1-f)
    X = X[, f >= maf.min, drop=FALSE]
  }
  keep = apply(X, 2, function(x) length(unique(x)) > 1)
  X = X[, keep, drop=FALSE]
  if (!ncol(X)) return(a)
  X = scale(X, center=TRUE, scale=FALSE)        # centre per SNP over all 808
  a$K = a$K + tcrossprod(X); a$m = a$m + ncol(X)
  a
}

t0 = Sys.time()
A_cand = acc(); for (g in CAND) A_cand = add(A_cand, g)
cat(sprintf("K_cand   %d SNPs  (%.1f min)\n", A_cand$m, as.numeric(difftime(Sys.time(),t0,units="mins")))); flush.console()

A_rest = acc(); A_maf = acc()
for (i in seq_along(rest)) {
  A_rest = add(A_rest, rest[i])
  A_maf  = add(A_maf,  rest[i], maf.min=0.01)
  if (i %% 1000 == 0) { cat(sprintf("  rest %5d/%d  SNPs=%d  (%.1f min)\n", i, length(rest),
      A_rest$m, as.numeric(difftime(Sys.time(),t0,units="mins")))); flush.console() }
}

mk = function(a) { K = a$K / a$m; dimnames(K) = list(ids, ids); K }
K_cand      = mk(A_cand)
K_rest      = mk(A_rest)
K_rest_maf01= mk(A_maf)
K_all       = { K = (A_cand$K + A_rest$K)/(A_cand$m + A_rest$m); dimnames(K)=list(ids,ids); K }

info = data.frame(kernel = c("K_cand","K_rest","K_rest_maf01","K_all"),
                  nSNP   = c(A_cand$m, A_rest$m, A_maf$m, A_cand$m+A_rest$m),
                  trace  = c(sum(diag(K_cand)), sum(diag(K_rest)), sum(diag(K_rest_maf01)), sum(diag(K_all))))
print(info)
save(K_cand, K_rest, K_rest_maf01, K_all, ids, info, file=file.path(OUT,"kernels.Rdata"))
write.csv(info, file.path(OUT,"kernel_info.csv"), row.names=FALSE)
cat(sprintf("\nDONE in %.1f min -> %s/kernels.Rdata\n", as.numeric(difftime(Sys.time(),t0,units="mins")), OUT))
