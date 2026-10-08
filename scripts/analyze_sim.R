# ---------------------------------------------------------------------------
# Analyse the revision simulation grid with the full method set.
#
# One SLURM array task per cell: simfun x nSNP x causal  (5 x 4 x 5 = 100).
# Each task evaluates all methods on all 100 replicates of that cell.
#
# Methods follow the original simulation1.R / simulation4.R exactly, plus BLUP
# (Reviewer 1, real-data comment 1).  LDpred2 and glinternet need packages that
# are not yet installed and are added in a second pass.
#
# Output: result/sim/eval_<simfun>_<nSNP>_<causal>.csv
#   columns: rep, method, mse, cor, secs
# (one CSV per cell rather than the original append-to-txt, so parallel array
#  tasks cannot interleave writes)
# ---------------------------------------------------------------------------
suppressMessages({library(MASS); library(glmnet); library(matrixcalc)})
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
for (f in list.files(paste0(RT,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
has.rrBLUP = requireNamespace("rrBLUP", quietly=TRUE)

nTrain=800; nTest=200; N=nTrain+nTest; nRegion=10
center_scale = function(x) scale(x, scale=FALSE)

SCEN = strsplit(Sys.getenv("SCENARIOS","linear,square,cosh,hyperbola,rickercurve"), "[,+]")[[1]]
grid = expand.grid(simfun=SCEN,
                   nSNP=c(20,100,500,1000), causal=c(0.2,0.4,0.6,0.8,1),
                   stringsAsFactors=FALSE)
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id>=1, id<=nrow(grid))
cf = grid[id,]; simfun=cf$simfun; nSNP=cf$nSNP; causal=cf$causal
REPS = as.integer(Sys.getenv("REPS","100"))
TAG  = Sys.getenv("DATA_TAG", "")

cell = sprintf("%sdata/sim%s/%s/%d/%s/", RT, TAG, simfun, nSNP, format(causal))
cat(sprintf("[task %d] %s nSNP=%d causal=%s  reps=%d\n", id, simfun, nSNP, causal, REPS))
flush.console()

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

for (k in seq_len(REPS)) {
  Xtr=list(); Xall=list()
  ok = TRUE
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

  # ---- our method: closed-form MINQUE -------------------------------------
  cfg = list(
    list("CLOSEFORM-MINQUE0-RIDGE", 5,   "MINQUE0", TRUE),
    list("CLOSEFORM-MINQUE1-RIDGE", 500, "MINQUE1", TRUE),
    list("CLOSEFORM-MINQUE0",       0,   "MINQUE0", FALSE),
    list("CLOSEFORM-MINQUE1",       0,   "MINQUE1", FALSE))
  for (cc in cfg) {
    r = tm(MINQUE.selection(ytr, Xtr, ik, MINQUE.type=cc[[3]], lambda=cc[[2]], constrain=cc[[4]]))
    p = if (r$ok) tryCatch(MINQUE.predict(ytr,Xall,ik,r$v$theta), error=function(e) NULL) else NULL
    add(k, cc[[1]], r, p, yts)
  }

  # ---- numerical MINQUE (ridge / lasso / elastic net) ----------------------
  for (a in list(c("NUMERICAL-MINQUE-RIDGE",0), c("NUMERICAL-MINQUE-LASSO",1),
                 c("NUMERICAL-MINQUE-ELASTICNET",0.5))) {
    r = tm(MINQUE.selection.numerical(ytr, Xtr, ik, alpha=as.numeric(a[2])))
    p = if (r$ok) tryCatch(MINQUE.predict(ytr,Xall,ik,r$v$theta), error=function(e) NULL) else NULL
    add(k, a[1], r, p, yts)
  }

  # ---- GMM-Lasso baseline --------------------------------------------------
  K.list = lapply(Xall, function(G) findKernel("product", G))
  yfull  = c(ytr, yts); idx = (nTrain+1):N
  r = tm(GMMLasso(y=yfull, index=idx, K=K.list))
  add(k, "GMM", r, if (r$ok) as.numeric(r$v$out[,1]) else NULL, yts)
  r = tm(GMMLasso(y=yfull, index=idx, K=K.list, lambda=0))
  add(k, "GMM-0", r, if (r$ok) as.numeric(r$v$out[,1]) else NULL, yts)

  # ---- linear BLUP (additive GRM over all regions) -------------------------
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
  if (k %% 5 == 0) { cat(sprintf("  rep %3d/%d\n", k, REPS)); flush.console() }
}

res = do.call(rbind, rows)
dir.create(paste0(RT,"result/sim",TAG), recursive=TRUE, showWarnings=FALSE)
write.csv(res, sprintf("%sresult/sim%s/eval_%s_%d_%s.csv", RT, TAG, simfun, nSNP, format(causal)), row.names=FALSE)
cat("[task ",id,"] DONE  rows=",nrow(res),"\n",sep="")
