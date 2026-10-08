# ---------------------------------------------------------------------------
# ADNI real-data analysis.  EXP=1 candidate genes only; EXP=2 adds the rest of
# the genome.  One array task per phenotype; all 100 splits inside it.
#
# Every method in both experiments reads the SAME stored split, so iteration k
# is the identical partition everywhere and the two experiments compare cell by
# cell rather than only on averages.
#
# Covariates (AGE, PTGENDER, PTEDUCAT) are fitted on TRAINING rows only each
# split; genetic models fit the training residual and the reported prediction is
# covariate + genetic, scored against the raw phenotype.  The published pipeline
# residualised on the whole sample, which leaked test information into training.
# ---------------------------------------------------------------------------
LIB="/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/Rlib"; .libPaths(c(LIB,.libPaths()))
suppressMessages({library(MASS); library(glmnet); library(Matrix)})
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
SH = paste0(RT,"realdata/shared/")
source(paste0(RT,"R/GMMLasso.R")); source(paste0(RT,"R/ChooseLambda.R")); source(paste0(RT,"R/ScaleKernel.R"))
has = function(p) requireNamespace(p, quietly=TRUE)

EXP  = as.integer(Sys.getenv("EXP","1"))
REPS = as.integer(Sys.getenv("REPS","100"))
BUDGET = as.numeric(Sys.getenv("BUDGET","900"))
ONLY   = Sys.getenv("ONLY","")          # "ldpred2" reruns only LDpred2, into <pheno>_ldpred2.csv
want.ld = ONLY %in% c("","ldpred2")

load(paste0(SH,"pheno_splits.Rdata"))                       # dat, splits, PHENO, COVAR
load(paste0(SH,"kernels/kernels.Rdata"))                    # K_cand,K_rest,K_rest_maf01,K_all,ids
load(paste0(SH,"kernels/kernels_candidate_individual.Rdata")) # K.cand.list, nsnp

id  = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); ph = PHENO[id]
OUT = paste0(RT,"realdata/", if (EXP==1) "exp1_candidate" else "exp2_genomewide", "/result/")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
cat(sprintf("[EXP%d] phenotype=%s reps=%d\n", EXP, ph, REPS)); flush.console()

# --- genotypes for LDpred2 -------------------------------------------------
G = NULL; LDSF = NULL
if (has("bigsnpr") && want.ld) {
  if (EXP==1) {
    CAND = names(K.cand.list); cols=list()
    for (g in CAND) { load(paste0(RT,"data/ADNI/",g,".Rdata")); cols[[g]] = as.matrix(geno) }
    G = do.call(cbind, cols); gene.of = rep(CAND, times=sapply(cols, ncol))
  } else {
    load(paste0(SH,"kernels/pruned_genomewide.Rdata")); G = G.pruned     # also supplies gene.of
  }
  sdg = apply(G,2,sd); G = G[, sdg>0, drop=FALSE]; gene.of = gene.of[sdg>0]
  # LD reference: block-diagonal, one block per gene, from all genotyped
  # individuals.  LD is a property of the genotypes, not the phenotype, so it
  # carries no test-set outcome information; only the summary statistics are
  # recomputed per split from training rows.  The first run formed a dense
  # cor() over every SNP -- 649,072^2 doubles, 3.1 TB -- and could not allocate.
  blk = split(seq_len(ncol(G)), factor(gene.of, levels=unique(gene.of)))
  ii = vector("list",length(blk)); jj = ii; xx = ii
  for (t in seq_along(blk)) { b = blk[[t]]
    C = if (length(b)>1) suppressWarnings(cor(G[,b,drop=FALSE])) else matrix(1)
    C[is.na(C)] = 0; C[abs(C) < 0.02] = 0; diag(C) = 1
    nz = which(C != 0, arr.ind=TRUE); ii[[t]] = b[nz[,1]]; jj[[t]] = b[nz[,2]]; xx[[t]] = C[nz] }
  LDmat = Matrix::sparseMatrix(i=unlist(ii), j=unlist(jj), x=unlist(xx), dims=c(ncol(G),ncol(G)))
  LDSF  = bigsparser::as_SFBM(methods::as(LDmat,"dgCMatrix"),
                              backingfile=tempfile("ld", tmpdir=Sys.getenv("TMPDIR")))
  rm(LDmat, ii, jj, xx); invisible(gc())
  cat(sprintf("LD reference: %d SNPs in %d gene blocks, %d non-zeros\n",
              ncol(G), length(blk), length(LDSF$p) )); flush.console()
}

# --- variance-component basis ---------------------------------------------
# interaction expansion over the supplied kernels, plus optional additive-only
# background terms that stay OUT of the expansion (the GRM idiom)
build.basis = function(K.list, background=list()) {
  n = nrow(K.list[[1]]); J = matrix(1,n,n)
  base = c(list(J), K.list)
  v = list(diag(1,n,n))
  for (i in seq_along(base)) for (j in i:length(base)) v[[length(v)+1]] = base[[i]]*base[[j]]
  c(v, background)
}
minque = function(y, v, lambda, type, constrain) {
  K = length(v); N = length(y)
  A = if (type=="MINQUE0") diag(1,N,N) else solve(Reduce(`+`, mapply(`*`, v, rep(1/K,K), SIMPLIFY=FALSE)))
  Rv = lapply(v, function(V) A %*% V); Ry = A %*% y
  u = vapply(v, function(V) sum(Ry * (V %*% Ry)), 0)
  Fm = matrix(0,K,K)
  for (i in 1:K) for (j in i:K) { Fm[i,j] = sum(Rv[[i]]*Rv[[j]]); Fm[j,i] = Fm[i,j] }
  pen = diag(lambda,K,K); pen[1,1] = 0
  th = as.numeric(MASS::ginv(Fm+pen) %*% u)
  if (constrain) th[th<0] = 0
  th
}
minque.pred = function(ytr, v, theta, ntr, nte) {
  KU = Reduce(`+`, mapply(`*`, v[-1], theta[-1], SIMPLIFY=FALSE))
  sig = KU[(ntr+1):(ntr+nte), 1:ntr]; Sig = KU[1:ntr,1:ntr] + diag(theta[1], ntr)
  as.numeric(sig %*% MASS::ginv(Sig) %*% ytr)
}
tm = function(e){ t0=Sys.time(); v=tryCatch({setTimeLimit(elapsed=BUDGET,transient=TRUE)
  on.exit(setTimeLimit(elapsed=Inf,transient=TRUE)); force(e)},
  error=function(er) structure(list(m=conditionMessage(er)),class="fx"))
  list(v=v,s=as.numeric(difftime(Sys.time(),t0,units="secs")),ok=!inherits(v,"fx"),
       msg=if(inherits(v,"fx")) v$m else NA_character_) }

rows=list()
add=function(k,meth,r,pred,ytest){ mse=NA_real_; cc=NA_real_
  if (r$ok && !is.null(pred) && all(is.finite(pred)) && sd(pred)>0){ mse=mean((ytest-pred)^2); cc=cor(ytest,pred) }
  rows[[length(rows)+1]] <<- data.frame(phenotype=ph,rep=k,method=meth,mse=mse,cor=cc,secs=r$s,msg=r$msg) }

for (k in seq_len(REPS)) {
  sp = splits[[ph]][[k]]; tr = sp$train; te = sp$test
  ntr=length(tr); nte=length(te); idx=c(tr,te)
  y = dat[[ph]]; X.cov = as.matrix(dat[, COVAR])

  # covariates fitted on training rows only
  cf  = lm.fit(cbind(1, X.cov[tr,,drop=FALSE]), y[tr])$coefficients; cf[is.na(cf)] = 0
  pc.tr = as.numeric(cbind(1,X.cov[tr,,drop=FALSE]) %*% cf)
  pc.te = as.numeric(cbind(1,X.cov[te,,drop=FALSE]) %*% cf)
  r.tr  = y[tr] - pc.tr; y.te = y[te]
  if (ONLY != "ldpred2") {
  add(k,"COVARIATE-ONLY",list(s=0,ok=TRUE,msg=NA),pc.te,y.te)

  Ksub = lapply(K.cand.list, function(K) scaleK(K[idx,idx]))
  bg   = if (EXP==2) list(scaleK(K_rest_maf01[idx,idx])) else list()
  v    = build.basis(Ksub, bg)
  v.tr = lapply(v, function(V) V[1:ntr,1:ntr])

  for (cfg in list(c("CF-MINQUE0-RIDGE",5,"MINQUE0",TRUE), c("CF-MINQUE1-RIDGE",500,"MINQUE1",TRUE),
                   c("CF-MINQUE0",0,"MINQUE0",FALSE),      c("CF-MINQUE1",0,"MINQUE1",FALSE))) {
    r = tm(minque(r.tr, v.tr, as.numeric(cfg[2]), cfg[3], as.logical(cfg[4])))
    p = if (r$ok) tryCatch(pc.te + minque.pred(r.tr, v, r$v, ntr, nte), error=function(e) NULL) else NULL
    add(k, cfg[1], r, p, y.te)
  }
  # sensitivity: background inside the interaction expansion instead of additive
  if (EXP==2) {
    v2 = build.basis(c(Ksub, list(scaleK(K_rest_maf01[idx,idx]))))
    r = tm(minque(r.tr, lapply(v2,function(V) V[1:ntr,1:ntr]), 500, "MINQUE1", TRUE))
    p = if (r$ok) tryCatch(pc.te + minque.pred(r.tr, v2, r$v, ntr, nte), error=function(e) NULL) else NULL
    add(k,"CF-MINQUE1-RIDGE-FULLINT", r, p, y.te)
  }
  # numerical MINQUE on the same basis
  for (a in list(c("NUM-RIDGE",0),c("NUM-LASSO",1),c("NUM-ENET",0.5))) {
    r = tm({ T = matrix(unlist(lapply(v.tr, c)), nrow=ntr*ntr, ncol=length(v.tr))
             pf = c(0, rep(1, ncol(T)-1))
             as.numeric(coef(cv.glmnet(T, c(r.tr%*%t(r.tr)), alpha=as.numeric(a[2]),
                                       lower.limits=0, penalty.factor=pf), s="lambda.min"))[-1] })
    p = if (r$ok) tryCatch(pc.te + minque.pred(r.tr, v, r$v, ntr, nte), error=function(e) NULL) else NULL
    add(k, a[1], r, p, y.te)
  }
  # GMM on the region kernels (plus background in EXP 2)
  Kg = c(Ksub, bg)
  r = tm(GMMLasso(y=c(r.tr, y.te-pc.te), index=(ntr+1):(ntr+nte), K=Kg))
  add(k,"GMM", r, if (r$ok) pc.te + as.numeric(r$v$out[,1]) else NULL, y.te)
  r = tm(GMMLasso(y=c(r.tr, y.te-pc.te), index=(ntr+1):(ntr+nte), K=Kg, lambda=0))
  add(k,"GMM-0", r, if (r$ok) pc.te + as.numeric(r$v$out[,1]) else NULL, y.te)

  # linear BLUP: candidate GRM in EXP 1, whole-genome GRM in EXP 2
  if (has("rrBLUP")) {
    Kb = scaleK((if (EXP==1) K_cand else K_all)[idx,idx])
    r = tm({ ms = rrBLUP::mixed.solve(y=r.tr, K=Kb[1:ntr,1:ntr]); lam = ms$Ve/ms$Vu
             as.numeric(Kb[(ntr+1):(ntr+nte),1:ntr] %*% solve(Kb[1:ntr,1:ntr]+lam*diag(ntr), r.tr)) })
    add(k,"BLUP", r, if (r$ok) pc.te + r$v else NULL, y.te)
  }
  }  # end non-LDpred2 methods
  # LDpred2
  if (!is.null(LDSF)) {
    r = tm({
      Gtr = G[tr,,drop=FALSE]; mu = colMeans(Gtr); sdv = apply(Gtr, 2, sd)
      okc = is.finite(sdv) & sdv > 0                   # polymorphic in THIS training set
      Xtr = sweep(sweep(Gtr,2,mu,`-`),2,ifelse(okc,sdv,1),`/`); rm(Gtr)
      sxx = colSums(Xtr^2); xty = as.numeric(crossprod(Xtr, r.tr))
      b   = ifelse(okc, xty/sxx, 0)
      se  = sqrt(pmax(sum(r.tr^2) - b*xty, 0)/(ntr-2)/sxx)
      # LDpred2 requires strictly positive SEs.  SNPs with no variation in this
      # training split enter as beta 0 with unit SE -- no signal -- which keeps
      # the LD reference dimension fixed across splits.  The first run passed
      # NaN SEs for exactly these SNPs and every fit was rejected.
      bad = !okc | !is.finite(se) | se <= 0; b[bad] = 0; se[bad] = 1
      df = data.frame(beta=b, beta_se=se, n_eff=ntr)
      bt = tryCatch({
             au = bigsnpr::snp_ldpred2_auto(LDSF, df, h2_init=0.2, vec_p_init=c(0.01,0.1),
                                            burn_in=100, num_iter=100, ncores=1)
             Bm = sapply(au, function(z) z$beta_est); keepc = apply(Bm, 2, function(z) all(is.finite(z)))
             if (!any(keepc)) stop("all chains diverged"); rowMeans(Bm[, keepc, drop=FALSE]) },
           error=function(e) bigsnpr::snp_ldpred2_inf(LDSF, df, h2=0.2))
      bt[!is.finite(bt)] = 0
      Xte = sweep(sweep(G[te,,drop=FALSE],2,mu,`-`),2,ifelse(okc,sdv,1),`/`); Xte[, !okc] = 0
      as.numeric(Xte %*% bt) })
    add(k,"LDpred2", r, if (r$ok) pc.te + r$v else NULL, y.te)
  }
  if (k %% 10 == 0) { cat(sprintf("  rep %3d/%d\n",k,REPS)); flush.console() }
}
res = do.call(rbind, rows)
write.csv(res, sprintf("%s%s%s.csv", OUT, ph, if (ONLY=="ldpred2") "_ldpred2" else ""), row.names=FALSE)
cat("[EXP",EXP,"] ",ph," DONE rows=",nrow(res),"\n",sep="")
