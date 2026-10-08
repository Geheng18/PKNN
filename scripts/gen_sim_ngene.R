# ---------------------------------------------------------------------------
# Data for the "how many candidate genes?" experiment (AE comment 1.4).
#
# Fixed:   h2 = 0.5, 20 SNPs per gene, 2 CAUSAL regions, N = 1000.
# Varied:  total number of input regions R = 2,4,6,8,10,12,16,20.
#
# Holding the causal count fixed at 2 while growing R is what answers the
# question that was asked.  R = 2 is the oracle -- only the informative genes
# are supplied -- and R = 20 supplies 18 uninformative ones alongside them, so
# the curve measures exactly the cost of over-inclusion.  Letting the causal
# count grow with R instead would confound dilution with a change in the
# amount of signal, which the existing causal-rate axis already covers.
#
# R = 20 is the ceiling: pick.regions places one region per chromosome, and
# bchQC03 contains 20 chromosomes (4 and 5 are absent), every one of which is
# large enough to hold a 20-SNP window.
#
# Output: data/sim_ngene/<simfun>/<R>/Sim_<rep>_<region>_<R>.Rdata , Phe_<rep>.txt
# ---------------------------------------------------------------------------
suppressMessages({library(BEDMatrix); library(data.table); library(MASS)})
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
for (f in list.files(paste0(RT,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
UKB = "/home/heng.ge/blue/heng.ge/KNN_Data_backup_again/data/UKB/bchQC03"

N          = 1000
nSNP       = 20
H2         = 0.5
nCausal    = 2
REGION.LIST = c(2,4,6,8,10,12,16,20)
simfun.list = c("linear","cosh","square","hyperbola","rickercurve")

grid = expand.grid(simfun = simfun.list, nRegion = REGION.LIST, stringsAsFactors = FALSE)
id   = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id >= 1, id <= nrow(grid))
cf   = grid[id,]; simfun = cf$simfun; nRegion = cf$nRegion
times = as.integer(Sys.getenv("REPS","100"))

cat(sprintf("[task %d] simfun=%s nRegion=%d causal=%d h2=%.2f reps=%d\n",
            id, simfun, nRegion, nCausal, H2, times)); flush.console()

out = file.path(RT, "data/sim_ngene", simfun, nRegion)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

bm  = BEDMatrix(paste0(UKB,".bed"), simple_names = TRUE)
bim = fread(paste0(UKB,".bim"), header = FALSE)

t0 = Sys.time()
for (k in seq_len(times)) {
  set.seed(20260930 + 1000*id + k)          # deterministic, unique per (cell, rep)
  d = Simdata.LD(bm, bim, N = N, nSNP = nSNP, nRegion = nRegion,
                 nRegion.true = nCausal, simfun = simfun, h2 = H2)
  for (l in seq_len(nRegion)) {
    geno = d$X[[l]]
    save(geno, file = file.path(out, sprintf("Sim_%d_%d_%d.Rdata", k, l, nRegion)))
  }
  write.table(d$y, file.path(out, sprintf("Phe_%d.txt", k)),
              row.names = FALSE, col.names = FALSE)
  write.csv(d$info, file.path(out, sprintf("regions_%d.csv", k)), row.names = FALSE)
  if (k %% 10 == 0) { cat(sprintf("  rep %3d/%d  (%.1f min)\n", k, times,
      as.numeric(difftime(Sys.time(), t0, units="mins")))); flush.console() }
}
cat("[task ", id, "] DONE\n", sep="")
