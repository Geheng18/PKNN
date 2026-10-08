# ---------------------------------------------------------------------------
# Convert ADNI gene-level PLINK filesets -> .Rdata matrices
#
# Source (READ-ONLY, never modified): /home/heng.ge/orange/lucienq/adni/QC/adni_gene
# Output:                             KNN_VC_Selection_Revision/data/ADNI/<gene>.Rdata
#
# Reproduces the original ADNI.R recipe exactly:
#   - genotypes read from .bed (additive coding 0/1/2)
#   - missing genotypes imputed by two random Bernoulli(MAF) allele draws
#   - rownames set to the PTID column of the .fam file
# Verified against the existing data/ADNI/*.Rdata: 100% agreement on all
# non-missing entries, identical dim / rownames / allele coding.
#
# plinkBED is not installed on HiPerGator; BEDMatrix 2.0.3 is used instead and
# gives byte-identical genotypes.
#
# Usage:  Rscript convert_adni.R [gene_list_file]
#         default gene list = the 11 candidate genes used in the paper
# ---------------------------------------------------------------------------
suppressMessages(library(BEDMatrix))

SRC  = "/home/heng.ge/orange/lucienq/adni/QC/adni_gene"
DEST = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/data/ADNI"
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

# impute missing genotypes from MAF (identical to `fill` in the original ADNI.R)
fill = function(SNP) {
  ind.na = which(is.na(SNP))
  if (length(ind.na) == 0) return(SNP)
  N = length(SNP); nN = length(ind.na)
  MAF = sum(SNP, na.rm = TRUE) / (2 * (N - nN))
  SNP[ind.na] = rbinom(nN, 1, MAF) + rbinom(nN, 1, MAF)
  SNP
}

candidate.genes = c("05_141971742","08_038854504","09_134378288","11_005274420",
                    "12_053901639","14_094843083","19_045417503","19_045409005",
                    "19_045394476","21_036160097","21_048018530")

args = commandArgs(trailingOnly = TRUE)
genes = if (length(args) >= 1) readLines(args[1]) else candidate.genes
genes = trimws(genes); genes = genes[nzchar(genes)]

set.seed(20260827)   # imputation is stochastic -- fixed for reproducibility

log = data.frame()
for (g in genes) {
  bed = file.path(SRC, paste0(g, ".bed"))
  if (!file.exists(bed)) { warning("missing fileset: ", g); next }

  bm  = BEDMatrix(bed, simple_names = TRUE)
  X   = as.matrix(bm)
  nNA = sum(is.na(X))

  fam = read.table(file.path(SRC, paste0(g, ".fam")), stringsAsFactors = FALSE)
  stopifnot(identical(rownames(X), as.character(fam$V2)))

  geno = apply(X, 2, fill)
  rownames(geno) = fam$V2
  stopifnot(!anyNA(geno))

  save(geno, file = file.path(DEST, paste0(g, ".Rdata")))
  log = rbind(log, data.frame(gene = g, n = nrow(geno), nSNP = ncol(geno),
                              pct.imputed = round(100 * nNA / length(X), 4)))
  cat(sprintf("%-16s %4d ind x %5d SNPs  (%.3f%% imputed)\n",
              g, nrow(geno), ncol(geno), 100 * nNA / length(X)))
}
write.csv(log, file.path(dirname(DEST), "ADNI_conversion_log.csv"), row.names = FALSE)
cat("\nconverted", nrow(log), "genes ->", DEST, "\n")
