# ---------------------------------------------------------------------------
# Fair lambda tuning for CF-MINQUE1-ridge on ADNI.
#
# Every other method in the real-data tables selects its tuning parameter from
# the training data (cv.glmnet for NUM and GMM, REML for BLUP, auto for LDpred2).
# CF-MINQUE1-ridge used a fixed lambda = 500, chosen for RAW kernels; the revision
# fed it scaleK'd kernels (x~18.5 per kernel, x~340 per interaction product), so
# that lambda became a far weaker penalty than the published one.
#
# Here lambda is chosen by 5-fold cross-validation INSIDE each training set --
# the test rows are never used -- over a wide log grid, with the tuning
# criterion fixed in advance as mean inner-fold MSE of the genetic residual.
# The fixed lambda = 500 arms are recorded from the same fits for reference.
# Both kernel scales are run: raw (as published) and scaled (as the revision).
# ---------------------------------------------------------------------------
suppressMessages({library(MASS)})
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"; SH = paste0(RT,"realdata/shared/")
source(paste0(RT,"R/ScaleKernel.R"))
load(paste0(SH,"pheno_splits.Rdata")); load(paste0(SH,"kernels/kernels_candidate_individual.Rdata"))
EXP = as.integer(Sys.getenv("EXP","1")); REPS = as.integer(Sys.getenv("REPS","100"))
if (EXP == 2) load(paste0(SH,"kernels/kernels.Rdata"))
PH4 = c("AV45","Entorhinal","FDG","Hippocampus")
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); ph = PH4[id]
OUT = paste0(RT,"realdata/lambda_cv/", if (EXP==1) "exp1_candidate" else "exp2_genomewide", "/result/")
GRID = c(0, 10^seq(-3, 9, by = 0.5))
NFOLD = 5
cat(sprintf("[lambda-cv EXP%d] %s reps=%d grid=%d values\n", EXP, ph, REPS, length(GRID))); flush.console()

build.basis = function(K.list) {
  n = nrow(K.list[[1]]); base = c(list(matrix(1,n,n)), K.list); v = list(diag(1,n,n))
  for (i in seq_along(base)) for (j in i:length(base)) v[[length(v)+1]] = base[[i]]*base[[j]]
  v }
# lambda-independent part of MINQUE1: A, R V_i, u, F  (computed once per fit)
prep = function(y, v) {
  K = length(v)
  A  = solve(Reduce(`+`, mapply(`*`, v, rep(1/K,K), SIMPLIFY=FALSE)))
  Rv = lapply(v, function(V) A %*% V); Ry = A %*% y
  u  = vapply(v, function(V) sum(Ry * (V %*% Ry)), 0)
  Fm = matrix(0,K,K); for (i in 1:K) for (j in i:K) { Fm[i,j] = sum(Rv[[i]]*Rv[[j]]); Fm[j,i] = Fm[i,j] }
  list(F=Fm, u=u, K=K) }
theta.at = function(p, lambda) { pen = diag(lambda,p$K,p$K); pen[1,1] = 0
  th = as.numeric(MASS::ginv(p$F + pen) %*% p$u); th[th<0] = 0; th }
pred = function(ytr, v, th, a, b) {          # a = fit rows, b = predict rows (indices into v)
  KU = Reduce(`+`, mapply(`*`, v[-1], th[-1], SIMPLIFY=FALSE))
  S  = KU[a,a] + diag(th[1], length(a))
  w  = tryCatch(solve(S, ytr), error=function(e) MASS::ginv(S) %*% ytr)
  as.numeric(KU[b,a] %*% w) }

rows = list()
for (k in seq_len(REPS)) {
  sp = splits[[ph]][[k]]; tr = sp$train; te = sp$test; ntr = length(tr); nte = length(te); idx = c(tr,te)
  y = dat[[ph]]; Xc = as.matrix(dat[, COVAR])
  cf = lm.fit(cbind(1, Xc[tr,,drop=FALSE]), y[tr])$coefficients; cf[is.na(cf)] = 0
  pc.tr = as.numeric(cbind(1,Xc[tr,,drop=FALSE]) %*% cf); pc.te = as.numeric(cbind(1,Xc[te,,drop=FALSE]) %*% cf)
  r.tr = y[tr] - pc.tr; r.te = y[te] - pc.te

  Kr = lapply(K.cand.list, function(K) K[idx,idx])
  if (EXP == 2) Kr = c(Kr, list(K_rest_maf01[idx,idx]))          # full interaction: K_rest inside the expansion
  for (scl in c("raw","scaled")) {
    Kin = if (scl == "raw") Kr else lapply(Kr, scaleK)
    v = build.basis(Kin)

    # inner 5-fold CV on the training rows only
    set.seed(1000*k + id)
    fold = sample(rep(seq_len(NFOLD), length.out = ntr)); cvmse = matrix(NA, NFOLD, length(GRID))
    for (f in seq_len(NFOLD)) {
      a = which(fold != f); b = which(fold == f)
      p = prep(r.tr[a], lapply(v, function(V) V[a,a]))
      for (g in seq_along(GRID)) { th = theta.at(p, GRID[g])
        cvmse[f,g] = mean((r.tr[b] - pred(r.tr[a], v, th, a, b))^2) }
    }
    lam = GRID[which.min(colMeans(cvmse))]

    # refit on the full training set, score on the untouched test rows
    p = prep(r.tr, lapply(v, function(V) V[1:ntr,1:ntr]))
    for (arm in list(c("CV", lam), c("fixed500", 500))) {
      th = theta.at(p, as.numeric(arm[2])); g = pred(r.tr, v, th, 1:ntr, (ntr+1):(ntr+nte))
      rows[[length(rows)+1]] = data.frame(phenotype=ph, rep=k, kernels=scl, lambda.rule=arm[1],
        lambda=as.numeric(arm[2]),
        cor_final   = cor(y[te], pc.te + g),
        cor_genetic = if (sd(g) > 0) cor(r.te, g) else NA_real_,
        collapsed   = sd(g) == 0)
    }
  }
  write.csv(do.call(rbind, rows), paste0(OUT, ph, ".csv"), row.names = FALSE)
  if (k %% 5 == 0) { cat(sprintf("  rep %3d/%d\n", k, REPS)); flush.console() }
}
cat("[lambda-cv] DONE\n")
