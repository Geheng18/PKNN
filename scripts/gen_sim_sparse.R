# ---------------------------------------------------------------------------
# Within-region sparsity simulation (Reviewer 2, comment 1).
#
# Fixed:   h2 = 0.5, 20 SNPs per gene, 10 regions, N = 1000.
# Varied:  n.causal.snp = 5, 10, 15, 20  -- how many SNPs inside each CAUSAL
#          region actually carry an effect
#          causal rate = 0.2 .. 1.0      -- how many of the 10 regions are causal
#          9 architectures               -- 5 link functions + 4 interaction
#
# The model always receives all 20 SNPs per region.  At n.causal.snp = 5, fifteen
# of the twenty are noise INSIDE a causal region, which is the sparse
# architecture the reviewer asked about; n.causal.snp = 20 reproduces the
# existing design and is the reference level.
#
# Output: data/sim_sparse/<scenario>/<ncausal>/<causal>/
# ---------------------------------------------------------------------------
suppressMessages({library(BEDMatrix); library(data.table); library(MASS)})
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
for (f in list.files(paste0(RT,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
UKB = "/home/heng.ge/blue/heng.ge/KNN_Data_backup_again/data/UKB/bchQC03"

N = 1000; nSNP = 20; nRegion = 10; H2 = 0.5
NCAUSAL = c(5,10,15,20)
NONLIN  = c("linear","cosh","square","hyperbola","rickercurve")
INTER   = c("within.region.threshold.TRUE","within.region.threshold.FALSE",
            "outside.region.threshold.TRUE","outside.region.threshold.FALSE")
SCEN    = c(NONLIN, INTER)
TRUE.LIST = 2:10                                    # causal regions -> rate /10

grid = expand.grid(scen = SCEN, ncausal = NCAUSAL,
                   nRegion.true = c(2,4,6,8,10), stringsAsFactors = FALSE)
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id>=1, id<=nrow(grid))
cf = grid[id,]; scen=cf$scen; ncausal=cf$ncausal; nRegion.true=cf$nRegion.true
causal = nRegion.true / nRegion
times  = as.integer(Sys.getenv("REPS","100"))

cat(sprintf("[task %d/%d] %s  n.causal.snp=%d  causal=%.1f  h2=%.2f  reps=%d\n",
            id, nrow(grid), scen, ncausal, causal, H2, times)); flush.console()

out = file.path(RT, "data/sim_sparse", scen, ncausal, format(causal))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
bm  = BEDMatrix(paste0(UKB,".bed"), simple_names = TRUE)
bim = fread(paste0(UKB,".bim"), header = FALSE)
is.inter = scen %in% INTER

t0 = Sys.time()
for (k in seq_len(times)) {
  set.seed(20261002 + 1000*id + k)
  d = if (is.inter)
        Simdata.LD.INT(bm, bim, N=N, nSNP=nSNP, nRegion=nRegion,
                       nRegion.true=nRegion.true, scenario=scen, h2=H2,
                       n.causal.snp=ncausal)
      else
        Simdata.LD(bm, bim, N=N, nSNP=nSNP, nRegion=nRegion,
                   nRegion.true=nRegion.true, simfun=scen, h2=H2,
                   n.causal.snp=ncausal)
  for (l in seq_len(nRegion)) {
    geno = d$X[[l]]
    save(geno, file = file.path(out, sprintf("Sim_%d_%d_%d.Rdata", k, l, nRegion)))
  }
  write.table(d$y, file.path(out, sprintf("Phe_%d.txt", k)),
              row.names=FALSE, col.names=FALSE)
  write.csv(d$info, file.path(out, sprintf("regions_%d.csv", k)), row.names=FALSE)
  if (k %% 20 == 0) { cat(sprintf("  rep %3d/%d (%.1f min)\n", k, times,
      as.numeric(difftime(Sys.time(),t0,units="mins")))); flush.console() }
}
cat("[task ", id, "] DONE\n", sep="")
