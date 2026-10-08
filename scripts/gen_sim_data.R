# ---------------------------------------------------------------------------
# Task 1 -- generate revision simulation data with TRUE LD and reduced effect size
#
# Grid (one SLURM array task per cell):
#   simfun : linear, square, cosh, hyperbola, rickercurve      (5)
#   nSNP   : 20, 100, 500, 1000                                (4)
#   causal : 0.2, 0.4, 0.6, 0.8, 1.0  (nRegion.true = 2..10)   (5)
#   = 100 cells x 100 replicates
#
# Genotypes: UK Biobank bchQC03, contiguous blocks, one region per chromosome.
# Output layout mirrors the original so downstream code is drop-in:
#   data/sim/<simfun>/<nSNP>/<causal>/Sim_<rep>_<region>_10.Rdata   (object `geno`)
#   data/sim/<simfun>/<nSNP>/<causal>/Phe_<rep>.txt
#   data/sim/<simfun>/<nSNP>/<causal>/regions_<rep>.csv
# ---------------------------------------------------------------------------
suppressMessages({library(BEDMatrix); library(MASS); library(data.table)})

RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
UKB = "/home/heng.ge/blue/heng.ge/KNN_Data_backup_again/data/UKB/bchQC03"
source(paste0(RT, "R/KernelPool.R"))
source(paste0(RT, "R/Simdata.LD.R"))

N        = 1000
nRegion  = 10
times    = 100
H2       = as.numeric(Sys.getenv("H2", "0.20"))   # heritability of the transformed signal
TAG      = Sys.getenv("DATA_TAG", "")             # "" = the original h2=0.2 grid

simfun.list = c("linear","square","cosh","hyperbola","rickercurve")
nSNP.list   = c(20,100,500,1000)
true.list   = c(2,4,6,8,10)

grid = expand.grid(simfun = simfun.list, nSNP = nSNP.list, nRegion.true = true.list,
                   stringsAsFactors = FALSE)

id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot(id >= 1, id <= nrow(grid))
cell = grid[id, ]
simfun = cell$simfun; nSNP = cell$nSNP; nRegion.true = cell$nRegion.true
causal = nRegion.true / nRegion

cat(sprintf("[task %d] simfun=%s nSNP=%d causal=%.1f h2=%.2f times=%d\n",
            id, simfun, nSNP, causal, H2, times)); flush.console()

out = file.path(RT, paste0("data/sim", TAG), simfun, nSNP, causal)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

bm  = BEDMatrix(paste0(UKB, ".bed"), simple_names = TRUE)
bim = fread(paste0(UKB, ".bim"), header = FALSE)

t0 = Sys.time()
for (k in seq_len(times)) {
  # deterministic, unique per (cell, replicate)
  set.seed(id * 100000L + k)
  d = Simdata.LD(bm, bim, N = N, nSNP = nSNP, nRegion = nRegion,
                 nRegion.true = nRegion.true, simfun = simfun, h2 = H2)

  for (l in seq_len(nRegion)) {
    geno = d$X[[l]]
    storage.mode(geno) = "integer"      # 0/1/2 -- compresses well
    save(geno, file = file.path(out, sprintf("Sim_%d_%d_%d.Rdata", k, l, nRegion)),
         compress = "gzip")
  }
  write.table(d$y, file = file.path(out, sprintf("Phe_%d.txt", k)),
              row.names = FALSE, col.names = FALSE)
  write.csv(d$info, file.path(out, sprintf("regions_%d.csv", k)), row.names = FALSE)

  if (k %% 10 == 0) {
    cat(sprintf("  rep %3d/%d  (%.1f s elapsed)\n", k, times,
                as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush.console()
  }
}
cat(sprintf("[task %d] DONE in %.1f min -> %s\n", id,
            as.numeric(difftime(Sys.time(), t0, units = "mins")), out))
