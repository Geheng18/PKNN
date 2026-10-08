# ---------------------------------------------------------------------------
# Real-data setup, shared by both experiments.
#   1. raw phenotypes + demographic covariates from ADNIMERGE (READ-ONLY on orange)
#   2. identical train/test splits, so every method sees the same data each iteration
#
# Phenotypes are kept RAW.  The published pipeline residualised y on
# AGE+PTGENDER+PTEDUCAT using the *whole* sample, which leaks test information
# into training and leaves a covariate-only baseline with nothing to predict.
# Here covariates are fitted per split on the training rows only.
# ---------------------------------------------------------------------------
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/realdata/shared/"
SRC = "/home/heng.ge/orange/lucienq/adni/raw/phe/ADNIMERGE.csv"
FAM = "/home/heng.ge/orange/lucienq/adni/QC/adni_gene/05_141971742.fam"
dir.create(RT, recursive=TRUE, showWarnings=FALSE)

phe = read.csv(SRC, stringsAsFactors=FALSE)
phe = phe[phe$VISCODE=="bl", ]
fam = read.table(FAM, stringsAsFactors=FALSE)
cat("ADNIMERGE baseline rows:", nrow(phe), " genotyped:", nrow(fam), "\n")

PHENO = c("FDG","AV45","Ventricles","Hippocampus","WholeBrain","Entorhinal")
COVAR = c("AGE","PTGENDER","PTEDUCAT")
keep  = c("PTID", COVAR, PHENO)
phe   = phe[, keep]
phe$PTGENDER = as.integer(factor(phe$PTGENDER))          # Female/Male -> 1/2

# align to genotype order; genotype rows are identical across all gene filesets
ord = match(fam$V2, phe$PTID)
dat = phe[ord, ]; rownames(dat) = fam$V2
dat$genotyped = !is.na(ord)
cat("matched to genotypes:", sum(dat$genotyped), "of", nrow(fam), "\n\n")

set.seed(20260914)
nTest = 200; TIMES = 100
splits = list(); summ = data.frame()
for (ph in PHENO) {
  ok = which(dat$genotyped & !is.na(dat[[ph]]) & complete.cases(dat[, COVAR]))
  N  = length(ok); nTrain = N - nTest
  s  = lapply(seq_len(TIMES), function(k) { p = sample(ok); list(train=p[1:nTrain], test=p[(nTrain+1):N]) })
  splits[[ph]] = s
  summ = rbind(summ, data.frame(phenotype=ph, N=N, nTrain=nTrain, nTest=nTest))
  cat(sprintf("%-12s N=%3d  nTrain=%3d  nTest=%d\n", ph, N, nTrain, nTest))
}
save(dat, splits, PHENO, COVAR, file=paste0(RT,"pheno_splits.Rdata"))
write.csv(summ, paste0(RT,"split_summary.csv"), row.names=FALSE)
cat("\nsaved ", RT, "pheno_splits.Rdata\n", sep="")
