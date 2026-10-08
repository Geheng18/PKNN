# ---------------------------------------------------------------------------
# Reproduce the PUBLISHED realdata.R calls verbatim, to find what made the
# published GMM score 0.07-0.12 on Entorhinal/FDG/Hippocampus.
#
# The first ablation showed GMMLasso's estimate/predict scale mismatch costs GMM
# only 0.01-0.03, so it cannot explain the published gap.  The remaining
# difference from that ablation is how genotypes are centred:
#   published : train and test blocks each centred on their OWN mean, then rbind
#   ablation 1: one global centring over all genotyped individuals
# Arm "published" runs the original R functions on the original genotype files
# exactly as realdata.R did; arm "globalcentre" changes only the centring.
# ---------------------------------------------------------------------------
suppressMessages({library(MASS); library(glmnet); library(matrixcalc)})
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"; SH = paste0(RT,"realdata/shared/")
OLD = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Simulation/data/ADNI/"      # read-only, the published inputs
for (f in c("KernelPool.R","KNN2LMM.R","MINQUE.selection.R","MINQUE.predict.R","MINQUE.selection.numerical.R",
            "MINQUE.selection.numerical.cv.R","GMMLasso.R","ChooseLambda.R","ScaleKernel.R")) source(paste0(RT,"R/",f))
load(paste0(SH,"pheno_splits.Rdata"))
PH4 = c("AV45","Entorhinal","FDG","Hippocampus")
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); ph = PH4[id]; REPS = as.integer(Sys.getenv("REPS","30"))
GENES = paste0(c("05_141971742","08_038854504","09_134378288","11_005274420","12_053901639","14_094843083",
                 "19_045417503","19_045409005","19_045394476","21_036160097","21_048018530"), ".Rdata")
geno.list = lapply(GENES, function(g) { load(paste0(OLD,g)); as.matrix(geno) })
po = read.table(paste0(OLD,"phe.",ph,".txt"), stringsAsFactors=FALSE)
m  = match(po$V1, rownames(dat)); y.old = rep(NA_real_, nrow(dat)); y.old[m[!is.na(m)]] = po$V2[!is.na(m)]
center_scale = function(x) scale(x, scale = FALSE)
nr = length(GENES); ik = rep(list(c("product")), nr)
sc = function(a,b) if (is.null(b) || !all(is.finite(b)) || sd(b)==0) NA_real_ else cor(a,b)
cat("[ablation2]", ph, "\n"); flush.console()

rows = list()
for (k in seq_len(REPS)) {
  sp = splits[[ph]][[k]]
  tr = sp$train[!is.na(y.old[sp$train])]; te = sp$test[!is.na(y.old[sp$test])]
  id.tr = rownames(dat)[tr]; id.te = rownames(dat)[te]
  ytr = scale(y.old[tr]); yte = scale(y.old[te]); y = c(ytr, yte); ntr = length(tr); nte = length(te)
  for (cen in c("published","globalcentre")) {
    Xtr = list(); Xall = list(); Kl = list()
    for (g in seq_len(nr)) {
      G = geno.list[[g]]
      if (cen == "published") {
        a = apply(G[id.tr,,drop=FALSE], 2, center_scale); b = apply(G[id.te,,drop=FALSE], 2, center_scale)
      } else {
        Gc = scale(G, center=TRUE, scale=FALSE); a = Gc[id.tr,,drop=FALSE]; b = Gc[id.te,,drop=FALSE]
      }
      Xtr[[g]] = a; X = rbind(a, b); Xall[[g]] = X; Kl[[g]] = findKernel("product", X)
    }
    t1 = tryCatch({ f = MINQUE.selection(ytr, Xtr, ik, lambda=500, MINQUE.type="MINQUE1", constrain=TRUE)
                    MINQUE.predict(ytr, Xall, ik, f$theta) }, error=function(e) NULL)
    t2 = tryCatch({ f = MINQUE.selection.numerical(ytr, Xtr, ik, alpha=1)
                    MINQUE.predict(ytr, Xall, ik, f$theta) }, error=function(e) NULL)
    t3 = tryCatch(as.numeric(GMMLasso(y=y, index=(ntr+1):(ntr+nte), K=Kl)$out[,1]), error=function(e) NULL)
    for (arm in list(list("CF-MINQUE1-ridge",t1), list("NUM-lasso",t2), list("GMM",t3)))
      rows[[length(rows)+1]] = data.frame(phenotype=ph, rep=k, centring=cen, method=arm[[1]],
                                          cor=sc(as.numeric(yte), arm[[2]]))
  }
  write.csv(do.call(rbind, rows), paste0(RT,"realdata/ablation2/result/",ph,".csv"), row.names=FALSE)
  cat(sprintf("  rep %d/%d\n", k, REPS)); flush.console()
}
cat("[ablation2] DONE\n")
