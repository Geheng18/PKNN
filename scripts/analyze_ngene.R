# ---------------------------------------------------------------------------
# "How many candidate genes?" experiment -- core methods (AE comment 1.4).
#
# Grid: 5 link functions x 8 region counts (2,4,6,8,10,12,16,20) = 40 tasks.
# Fixed throughout: h2 = 0.5, 20 SNPs per gene, 2 causal regions, N = 1000.
#
# Method set is the five that appear in the reported figures.  The variants
# that are computed by analyze_sim.R but never plotted (CLOSEFORM-MINQUE0*,
# the unconstrained closed forms, NUMERICAL-MINQUE-RIDGE, GMM-0) are omitted:
# the interaction basis reaches 232 components at R = 20, where one closed-form
# fit costs 109 s, so carrying the unused variants would multiply the grid by
# about 2.5 for nothing that is reported.
#
# Output: result/sim_ngene/eval_<simfun>_<R>.csv
# ---------------------------------------------------------------------------
suppressMessages({library(MASS); library(glmnet); library(matrixcalc)})
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
for (f in list.files(paste0(RT,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
has.rrBLUP = requireNamespace("rrBLUP", quietly=TRUE)

nTrain=800; nTest=200; N=nTrain+nTest
center_scale = function(x) scale(x, scale=FALSE)

SCEN   = strsplit(Sys.getenv("SCENARIOS","linear+square+cosh+hyperbola+rickercurve"), "[,+]")[[1]]
REGION = as.integer(strsplit(Sys.getenv("REGIONS","2+4+6+8+10+12+16+20"), "[,+]")[[1]])
grid = expand.grid(simfun=SCEN, nRegion=REGION, stringsAsFactors=FALSE)
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id>=1, id<=nrow(grid))
cf = grid[id,]; simfun=cf$simfun; nRegion=cf$nRegion
REPS = as.integer(Sys.getenv("REPS","100"))

cell = sprintf("%sdata/sim_ngene/%s/%d/", RT, simfun, nRegion)
nComp = 2 + nRegion + nRegion*(nRegion+1)/2
cat(sprintf("[task %d] %s nRegion=%d (%.0f components) reps=%d\n",
            id, simfun, nRegion, nComp, REPS)); flush.console()

tm = function(expr) { t0=Sys.time()
  v = tryCatch(force(expr), error=function(e) structure(list(m=conditionMessage(e)),class="fx"))
  list(v=v, s=as.numeric(difftime(Sys.time(),t0,units="secs")), ok=!inherits(v,"fx")) }

rows = list()
add = function(rep, method, r, pred, yts) {
  mse = NA_real_; cc = NA_real_
  if (r$ok && !is.null(pred) && all(is.finite(pred)) && sd(pred) > 0) {
    mse = mean((yts-pred)^2); cc = cor(yts, pred)
  }
  rows[[length(rows)+1]] <<- data.frame(rep=rep, method=method, mse=mse, cor=cc, secs=r$s)
}
flush.rows = function() {
  if (!length(rows)) return(invisible())
  dir.create(paste0(RT,"result/sim_ngene"), recursive=TRUE, showWarnings=FALSE)
  write.csv(do.call(rbind, rows),
            sprintf("%sresult/sim_ngene/eval_%s_%d.csv", RT, simfun, nRegion), row.names=FALSE)
}

for (k in seq_len(REPS)) {
  Xtr=list(); Xall=list(); ok=TRUE
  for (l in 1:nRegion) {
    fn = sprintf("%sSim_%d_%d_%d.Rdata", cell, k, l, nRegion)
    if (!file.exists(fn)) { ok=FALSE; break }
    load(fn); G = as.matrix(geno)
    Xtr[[l]]  = apply(G[1:nTrain,],2,center_scale)
    Xall[[l]] = rbind(apply(G[1:nTrain,],2,center_scale), apply(G[(nTrain+1):N,],2,center_scale))
  }
  if (!ok) { cat("  rep",k,"missing files - skipped\n"); next }
  y = unlist(read.table(sprintf("%sPhe_%d.txt", cell, k)))
  ytr = scale(y[1:nTrain]); yts = scale(y[(nTrain+1):N])
  ik  = rep(list(c("product")), nRegion)

  # ---- proposed: closed-form MINQUE, ridge --------------------------------
  r = tm(MINQUE.selection(ytr, Xtr, ik, MINQUE.type="MINQUE1", lambda=500, constrain=TRUE))
  p = if (r$ok) tryCatch(MINQUE.predict(ytr,Xall,ik,r$v$theta), error=function(e) NULL) else NULL
  add(k, "CLOSEFORM-MINQUE1-RIDGE", r, p, yts)

  # ---- numerical MINQUE (lasso / elastic net) -----------------------------
  for (a in list(c("NUMERICAL-MINQUE-LASSO",1), c("NUMERICAL-MINQUE-ELASTICNET",0.5))) {
    r = tm(MINQUE.selection.numerical(ytr, Xtr, ik, alpha=as.numeric(a[2])))
    p = if (r$ok) tryCatch(MINQUE.predict(ytr,Xall,ik,r$v$theta), error=function(e) NULL) else NULL
    add(k, a[1], r, p, yts)
  }

  # ---- GMM-Lasso baseline -------------------------------------------------
  K.list = lapply(Xall, function(G) findKernel("product", G))
  yfull  = c(ytr, yts); idx = (nTrain+1):N
  r = tm(GMMLasso(y=yfull, index=idx, K=K.list))
  add(k, "GMM", r, if (r$ok) as.numeric(r$v$out[,1]) else NULL, yts)

  # ---- linear BLUP (additive GRM over all supplied regions) ---------------
  if (has.rrBLUP) {
    r = tm({
      Kall = Reduce(`+`, K.list); Kall = scaleK(Kall)
      Ktt = Kall[1:nTrain,1:nTrain]; Kts = Kall[idx,1:nTrain]
      ms  = rrBLUP::mixed.solve(y=as.numeric(ytr), K=Ktt)
      lam = ms$Ve/ms$Vu
      as.numeric(Kts %*% solve(Ktt + lam*diag(nTrain), as.numeric(ytr)))
    })
    add(k, "BLUP", r, if (r$ok) r$v else NULL, yts)
  }
  rm(K.list); gc(full=TRUE)
  # write as we go: a task killed near the wall clock keeps its finished reps
  if (k %% 5 == 0) { flush.rows(); cat(sprintf("  rep %3d/%d\n", k, REPS)); flush.console() }
}
flush.rows()
cat("[task ",id,"] DONE  rows=",length(rows),"\n",sep="")
