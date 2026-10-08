# Generate the four interaction scenarios into the existing simulation trees.
# 4 scenarios x 4 SNP densities x 5 causal rates x 2 heritability arms = 160 cells.
suppressMessages({library(BEDMatrix); library(MASS); library(data.table)})
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
UKB = "/home/heng.ge/blue/heng.ge/KNN_Data_backup_again/data/UKB/bchQC03"
source(paste0(RT,"R/KernelPool.R")); source(paste0(RT,"R/Simdata.LD.R")); source(paste0(RT,"R/Simdata.LD.INT.R"))
N = 1000; nRegion = 10; times = 100
SCEN = c("within.region.threshold.TRUE","within.region.threshold.FALSE",
         "outside.region.threshold.TRUE","outside.region.threshold.FALSE")
grid = do.call(rbind, lapply(list(c("","0.20"), c("_h2_0.5","0.50")), function(a)
  cbind(expand.grid(scenario=SCEN, nSNP=c(20,100,500,1000), nRegion.true=c(2,4,6,8,10),
                    stringsAsFactors=FALSE), TAG=a[1], H2=as.numeric(a[2]))))
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id>=1, id<=nrow(grid))
cf = grid[id,]; causal = cf$nRegion.true/nRegion
cat(sprintf("[%d] %s nSNP=%d causal=%.1f h2=%.2f arm=%s\n", id, cf$scenario, cf$nSNP, causal, cf$H2,
            ifelse(cf$TAG=="","0.2","0.5"))); flush.console()
out = file.path(RT, paste0("data/sim", cf$TAG), cf$scenario, cf$nSNP, format(causal))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
bm = BEDMatrix(paste0(UKB,".bed"), simple_names=TRUE); bim = fread(paste0(UKB,".bim"), header=FALSE)
t0 = Sys.time()
for (k in seq_len(times)) {
  set.seed(500000L + id*1000L + k)
  d = Simdata.LD.INT(bm, bim, N=N, nSNP=cf$nSNP, nRegion=nRegion,
                     nRegion.true=cf$nRegion.true, scenario=cf$scenario, h2=cf$H2)
  for (l in seq_len(nRegion)) { geno = d$X[[l]]; storage.mode(geno) = "integer"
    save(geno, file=file.path(out, sprintf("Sim_%d_%d_%d.Rdata", k, l, nRegion)), compress="gzip") }
  write.table(d$y, file.path(out, sprintf("Phe_%d.txt", k)), row.names=FALSE, col.names=FALSE)
  write.csv(d$info, file.path(out, sprintf("regions_%d.csv", k)), row.names=FALSE)
  if (k %% 20 == 0) { cat(sprintf("  rep %3d/%d (%.1f min)\n", k, times,
      as.numeric(difftime(Sys.time(),t0,units="mins")))); flush.console() }
}
cat(sprintf("[%d] DONE %.1f min -> %s\n", id, as.numeric(difftime(Sys.time(),t0,units="mins")), out))
