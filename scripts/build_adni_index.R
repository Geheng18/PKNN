# Index the 24,219 ADNI gene filesets: name, SNP count, physical span.
# Also measures within-gene LD by contiguous block size, for the Methods table.
suppressMessages({library(BEDMatrix); library(data.table)})
D = "/home/heng.ge/orange/lucienq/adni/QC/adni_gene"
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
beds = list.files(D, pattern = "\\.bed$", full.names = FALSE)
genes = sub("\\.bed$", "", beds)
cat("indexing", length(genes), "genes\n"); flush.console()
# nSNP from .bed size: 3 + ceil(N/4)*nSNP bytes, N = 808 -> 202 bytes/SNP
sz = file.size(file.path(D, beds))
idx = data.frame(gene = genes, nSNP = as.integer((sz - 3) / 202), stringsAsFactors = FALSE)
idx = idx[idx$nSNP > 0, ]
write.csv(idx, paste0(RT, "data/adni_gene_index.csv"), row.names = FALSE)
cat("wrote index:", nrow(idx), "genes\n")
for (t in c(20,100,500,1000)) cat(sprintf("  genes with >= %5d SNPs : %6d\n", t, sum(idx$nSNP >= t)))

# LD by contiguous block size within real ADNI genes
set.seed(1); cat("\nADNI within-gene LD (contiguous blocks, 808 individuals)\n")
cat(sprintf("%-8s %-12s %-14s %s\n", "nSNP", "mean r^2", "frac r^2>0.1", "span(kb)"))
bim.span = function(g, s, K) { b = fread(file.path(D, paste0(g,".bim")), header=FALSE)
                               diff(range(b$V4[s:(s+K-1)]))/1000 }
for (K in c(20,100,500,1000)) {
  r2s = c(); sp = c()
  cand = idx$gene[idx$nSNP >= K]
  for (b in 1:10) {
    g = sample(cand, 1)
    bm = BEDMatrix(file.path(D, paste0(g,".bed")), simple_names=TRUE)
    s = sample(ncol(bm)-K+1, 1)
    X = bm[, s:(s+K-1)]
    X = X[, apply(X,2,function(x) length(unique(x[!is.na(x)]))>1), drop=FALSE]
    r = suppressWarnings(cor(X, use="pairwise.complete.obs"))
    r2s = c(r2s, r[upper.tri(r)]^2); sp = c(sp, bim.span(g,s,K))
  }
  cat(sprintf("%-8d %-12.4f %-14.4f %.1f\n", K, mean(r2s,na.rm=TRUE), mean(r2s>0.1,na.rm=TRUE), median(sp)))
}
