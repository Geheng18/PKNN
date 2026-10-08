# ---------------------------------------------------------------------------
# Why do the revision ADNI results disagree with the published ones?
#
# Three things changed between pipelines, and this job varies them one at a time
# on the same stored splits for the four reported phenotypes:
#
#  1. Kernel scale.  Raw product kernels have mean diagonal ~0.054, so scaleK
#     multiplies them by ~18.5 (Hadamard products by ~340).
#       - CF-MINQUE: the published run used raw kernels; the revision used scaled
#         kernels, which changes what lambda = 500 means by orders of magnitude.
#       - GMMLasso estimates on scaleK'd kernels but predicts with the kernels it
#         was given.  The published run gave it raw kernels -> estimate/predict
#         mismatch.  The revision gave it scaled kernels -> consistent.
#       - Numerical MINQUE scaleK's each component to estimate, then
#         MINQUE.predict rebuilds raw components -> the same mismatch.
#  2. Phenotype.  Published: residualised on AGE/SEX/EDUC using the full sample.
#     Revision: raw phenotype, covariates fitted on training rows only.
#  3. Numerical lambda rule.  Published: MINQUE.selection.numerical.cv (1-SD
#     band); revision: lambda.min.
#
# Arms whose label ends in _old reproduce the published convention, _new the
# revision's; _fixonly changes only the kernel scale.
# ---------------------------------------------------------------------------
LIB="/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/Rlib"; .libPaths(c(LIB,.libPaths()))
suppressMessages({library(MASS); library(glmnet)})
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"; SH = paste0(RT,"realdata/shared/")
OLD = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Simulation/data/ADNI/"     # read-only
for (f in c("GMMLasso.R","ChooseLambda.R","ScaleKernel.R","MINQUE.selection.numerical.cv.R")) source(paste0(RT,"R/",f))
load(paste0(SH,"pheno_splits.Rdata")); load(paste0(SH,"kernels/kernels_candidate_individual.Rdata"))

PH4 = c("AV45","Entorhinal","FDG","Hippocampus")
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); ph = PH4[id]
REPS = as.integer(Sys.getenv("REPS","30"))
OUT = paste0(RT,"realdata/ablation/result/")
cat("[ablation]", ph, "reps", REPS, "\n"); flush.console()

po = read.table(paste0(OLD,"phe.",ph,".txt"), stringsAsFactors=FALSE)   # published residualised phenotype
m  = match(po$V1, rownames(dat)); y.old = rep(NA_real_, nrow(dat)); y.old[m[!is.na(m)]] = po$V2[!is.na(m)]

build.basis = function(K.list) {
  n = nrow(K.list[[1]]); base = c(list(matrix(1,n,n)), K.list); v = list(diag(1,n,n))
  for (i in seq_along(base)) for (j in i:length(base)) v[[length(v)+1]] = base[[i]]*base[[j]]
  v }
minque = function(y, v, lambda, type, constrain) {
  K = length(v); N = length(y)
  A = if (type=="MINQUE0") diag(1,N,N) else solve(Reduce(`+`, mapply(`*`, v, rep(1/K,K), SIMPLIFY=FALSE)))
  Rv = lapply(v, function(V) A %*% V); Ry = A %*% y
  u = vapply(v, function(V) sum(Ry * (V %*% Ry)), 0)
  Fm = matrix(0,K,K); for (i in 1:K) for (j in i:K) { Fm[i,j] = sum(Rv[[i]]*Rv[[j]]); Fm[j,i] = Fm[i,j] }
  pen = diag(lambda,K,K); pen[1,1] = 0
  th = as.numeric(MASS::ginv(Fm+pen) %*% u); if (constrain) th[th<0] = 0; th }
minque.pred = function(ytr, v, theta, ntr, nte) {
  KU = Reduce(`+`, mapply(`*`, v[-1], theta[-1], SIMPLIFY=FALSE))
  as.numeric(KU[(ntr+1):(ntr+nte),1:ntr] %*% MASS::ginv(KU[1:ntr,1:ntr] + diag(theta[1],ntr)) %*% ytr) }
num.theta = function(vtr, y, rule) {
  T  = matrix(unlist(lapply(vtr, c)), nrow=length(y)^2, ncol=length(vtr)); yy = c(y %*% t(y))
  pf = c(0, rep(1, ncol(T)-1))
  cvf = cv.glmnet(T, yy, alpha=1, lower.limits=0, penalty.factor=pf)
  if (rule=="old") { lam = MINQUE.selection.numerical.cv(cvf, 1)
                     as.numeric(glmnet(T, yy, alpha=1, lambda=lam, lower.limits=0, penalty.factor=pf)$beta) }
  else as.numeric(coef(cvf, s="lambda.min"))[-1] }
safe = function(e) tryCatch(e, error=function(er) NULL)

# genetic predictions on the test rows for every arm, given training rows first
arms = function(tr, te, ytr, ytest.for.gmm) {
  idx = c(tr,te); ntr = length(tr); nte = length(te)
  Kr = lapply(K.cand.list, function(K) K[idx,idx])                 # raw (published)
  Ks = lapply(Kr, scaleK)                                          # scaled (revision)
  vr = build.basis(Kr); vs = build.basis(Ks)
  vr.tr = lapply(vr, function(V) V[1:ntr,1:ntr]); vs.tr = lapply(vs, function(V) V[1:ntr,1:ntr])
  P = list()
  P$CF_old  = safe(minque.pred(ytr, vr, minque(ytr, vr.tr, 500, "MINQUE1", TRUE), ntr, nte))
  P$CF_new  = safe(minque.pred(ytr, vs, minque(ytr, vs.tr, 500, "MINQUE1", TRUE), ntr, nte))
  P$NUM_old     = safe(minque.pred(ytr, vr, num.theta(lapply(vr.tr, scaleK), ytr, "old"), ntr, nte))  # estimate scaled, predict raw
  P$NUM_fixonly = safe(minque.pred(ytr, vs, num.theta(vs.tr, ytr, "old"), ntr, nte))
  P$NUM_new     = safe(minque.pred(ytr, vs, num.theta(vs.tr, ytr, "min"), ntr, nte))
  yg = c(ytr, ytest.for.gmm); ix = (ntr+1):(ntr+nte)
  g  = safe(GMMLasso(y=yg, index=ix, K=Kr)); P$GMM_old = if (is.null(g)) NULL else as.numeric(g$out[,1])   # raw in -> mismatch
  Kc = lapply(Kr, function(K) K * ntr/sum(diag(K[1:ntr,1:ntr])))   # internal scaleK becomes a no-op
  g  = safe(GMMLasso(y=yg, index=ix, K=Kc)); P$GMM_new = if (is.null(g)) NULL else as.numeric(g$out[,1])
  P }
sc = function(a,b) if (is.null(b) || !all(is.finite(b)) || sd(b)==0) NA_real_ else cor(a,b)

rows = list()
for (k in seq_len(REPS)) {
  sp = splits[[ph]][[k]]; tr = sp$train; te = sp$test

  # (A) published phenotype convention: full-sample residuals, each block scaled separately
  tro = tr[!is.na(y.old[tr])]; teo = te[!is.na(y.old[te])]
  ya = as.numeric(scale(y.old[tro])); yb = as.numeric(scale(y.old[teo]))
  PA = arms(tro, teo, ya, yb)
  for (nm in names(PA)) rows[[length(rows)+1]] = data.frame(phenotype=ph, rep=k, pheno="published-residual",
                          arm=nm, metric="cor", value=sc(yb, PA[[nm]]))

  # (B) revision convention: raw phenotype, covariates on training rows only
  y = dat[[ph]]; Xc = as.matrix(dat[, COVAR])
  cf = lm.fit(cbind(1, Xc[tr,,drop=FALSE]), y[tr])$coefficients; cf[is.na(cf)] = 0
  pc.tr = as.numeric(cbind(1,Xc[tr,,drop=FALSE]) %*% cf); pc.te = as.numeric(cbind(1,Xc[te,,drop=FALSE]) %*% cf)
  r.tr = y[tr] - pc.tr; r.te = y[te] - pc.te
  PB = arms(tr, te, r.tr, r.te)
  rows[[length(rows)+1]] = data.frame(phenotype=ph, rep=k, pheno="revision-raw+cov", arm="COVARIATE-ONLY", metric="cor_final", value=cor(y[te], pc.te))
  for (nm in names(PB)) {
    g = PB[[nm]]
    rows[[length(rows)+1]] = data.frame(phenotype=ph, rep=k, pheno="revision-raw+cov", arm=nm, metric="cor_genetic", value=sc(r.te, g))
    rows[[length(rows)+1]] = data.frame(phenotype=ph, rep=k, pheno="revision-raw+cov", arm=nm, metric="cor_final",
                                        value=if (is.null(g)) NA_real_ else cor(y[te], pc.te + g))
  }
  write.csv(do.call(rbind, rows), paste0(OUT, ph, ".csv"), row.names=FALSE)
  cat(sprintf("  rep %d/%d\n", k, REPS)); flush.console()
}
cat("[ablation]", ph, "DONE\n")
