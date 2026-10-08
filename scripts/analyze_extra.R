# ---------------------------------------------------------------------------
# Reviewer-requested comparators, run over the existing simulation grids.
# Current results are untouched: output goes to result/sim<TAG>_extra/.
#
#   LDpred2     (bigsnpr)     - Reviewer 2, comment 4: Bayesian linear PRS
#   glinternet  (glinternet)  - Reviewer 1, sim comment 5: explicit pairwise
#                               interactions under strong hierarchy
#   hierNet     (hierNet)     - same family; maintained stand-in for SHIM,
#                               which was never on CRAN (GitHub only, 2017)
#   shim        (GitHub)      - attempted; used if the install succeeded
#
# Reviewer 1 asked for these "if computationally feasible".  Interaction
# searches are O(p^2), so each method is given a wall-clock budget per fit and
# records why it stopped.  Those refusals are a reportable result, not a gap.
#
# One array task per (arm, cell): 2 x 100 = 200.
# ---------------------------------------------------------------------------
LIB="/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/Rlib"; .libPaths(c(LIB,.libPaths()))
suppressMessages({library(MASS); library(Matrix)})
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
has = function(p) requireNamespace(p, quietly=TRUE)
BUDGET = as.numeric(Sys.getenv("BUDGET","600"))   # seconds per fit
# Dimension ceilings, set from the smoke tests.  Both interaction searches are
# O(p^2) in the number of predictors p = nRegion * nSNP: glinternet OOM'd at
# 64 GB with p = 10,000, and hierNet needed >7 min for 2 replicates at p = 200.
# Above these ceilings the method is recorded as infeasible rather than run --
# which is the direct answer to Reviewer 1's "if computationally feasible".
P_MAX_GLINTERNET = as.numeric(Sys.getenv("P_MAX_GLINTERNET","2000"))
P_MAX_HIERNET    = as.numeric(Sys.getenv("P_MAX_HIERNET","500"))
# glinternet tuning.  screenLimit caps how many main effects may spawn
# interaction candidates; it was set to 50 (of 200) purely for speed, and the
# degenerate-fit rate rises with the causal rate (13% -> 72%), so the screen is
# a prime suspect.  Exposed here so it can be tested rather than assumed.
SCREENLIMIT = Sys.getenv("SCREENLIMIT","50")
SCREENLIMIT = if (SCREENLIMIT %in% c("none","NA","")) NULL else as.numeric(SCREENLIMIT)
NFOLDS      = as.integer(Sys.getenv("NFOLDS","2"))
# Which comparators this task runs.  LDpred2 is cheap and runs over the whole
# grid; the two interaction searches are expensive enough to need their own
# submission on the cells where they can finish at all.
WANT = strsplit(Sys.getenv("METHODS","ldpred2+glinternet+hiernet"), "[,+]")[[1]]
want = function(m) m %in% WANT

nTrain=800; nTest=200; N=nTrain+nTest; nRegion=10
center_scale = function(x) scale(x, scale=FALSE)

arms = c("", "_h2_0.5")
SCEN = strsplit(Sys.getenv("SCENARIOS","linear,square,cosh,hyperbola,rickercurve"), "[,+]")[[1]]
grid = expand.grid(simfun=SCEN,
                   nSNP=c(20,100,500,1000), causal=c(0.2,0.4,0.6,0.8,1),
                   stringsAsFactors=FALSE)
full = do.call(rbind, lapply(arms, function(a) cbind(grid, TAG=a, stringsAsFactors=FALSE)))

id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id>=1, id<=nrow(full))
cf = full[id,]; simfun=cf$simfun; nSNP=cf$nSNP; causal=cf$causal; TAG=cf$TAG
REPS = as.integer(Sys.getenv("REPS","100"))
cell = sprintf("%sdata/sim%s/%s/%d/%s/", RT, TAG, simfun, nSNP, format(causal))
cat(sprintf("[%d] arm=%s %s P=%d causal=%s reps=%d\n", id,
            ifelse(TAG=="","h2=0.2","h2=0.5"), simfun, nSNP, causal, REPS)); flush.console()

# run an expression under a wall-clock budget; classify the outcome
budgeted = function(expr) {
  t0=Sys.time()
  v = tryCatch({ setTimeLimit(elapsed=BUDGET, transient=TRUE)
                 on.exit(setTimeLimit(elapsed=Inf, transient=TRUE))
                 force(expr) },
      error=function(e){ m=conditionMessage(e)
        structure(list(m=if(grepl("reached elapsed time limit|reached CPU time limit",m))
                          "over time budget" else m), class="fx") })
  list(v=v, s=as.numeric(difftime(Sys.time(),t0,units="secs")), ok=!inherits(v,"fx"),
       msg=if(inherits(v,"fx")) v$m else NA_character_)
}

rows=list()
add=function(rep,method,r,pred,yts,note=NA){
  mse=NA_real_; cc=NA_real_
  if (r$ok && !is.null(pred) && all(is.finite(pred)) && stats::sd(pred)>0) {
    mse=mean((yts-pred)^2); cc=stats::cor(yts,pred) }
  rows[[length(rows)+1]] <<- data.frame(rep=rep, method=method, mse=mse, cor=cc,
                                        secs=r$s, ok=r$ok,
                                        msg=ifelse(is.na(note), r$msg, note))
}

for (k in seq_len(REPS)) {
  Xtr=list(); Xte=list(); okf=TRUE
  for (l in 1:nRegion) {
    fn=sprintf("%sSim_%d_%d_%d.Rdata",cell,k,l,nRegion)
    if(!file.exists(fn)){okf=FALSE;break}
    load(fn); G=as.matrix(geno)
    Xtr[[l]]=apply(G[1:nTrain,],2,center_scale); Xte[[l]]=apply(G[(nTrain+1):N,],2,center_scale)
  }
  if(!okf){cat("  rep",k,"missing - skip\n");next}
  Atr=do.call(cbind,Xtr); Ate=do.call(cbind,Xte); p=ncol(Atr)
  y=unlist(read.table(sprintf("%sPhe_%d.txt",cell,k)))
  ytr=as.numeric(scale(y[1:nTrain])); yts=as.numeric(scale(y[(nTrain+1):N]))
  sdv=apply(Atr,2,stats::sd); keep=which(sdv>0)
  Atr=Atr[,keep,drop=FALSE]; Ate=Ate[,keep,drop=FALSE]; p=ncol(Atr)

  # ---- LDpred2 ------------------------------------------------------------
  if (has("bigsnpr") && want("ldpred2")) {
    r = budgeted({
      n=nTrain
      b=as.numeric(crossprod(Atr,ytr)/colSums(Atr^2))                 # marginal betas
      res=ytr - sweep(Atr,2,b,`*`)                                     # per-SNP residual scale
      se=sqrt(colSums(res^2)/(n-2)/colSums(Atr^2))
      df=data.frame(beta=b, beta_se=se, n_eff=n)
      # block-diagonal LD: one block per region (regions sit on distinct chromosomes)
      off=0; blocks=list()
      for (l in 1:nRegion) {
        m=ncol(Xtr[[l]]); idxl=intersect(keep, (off+1):(off+m))-off
        if(length(idxl)>1) blocks[[length(blocks)+1]]=stats::cor(Xtr[[l]][,idxl,drop=FALSE])
        off=off+m
      }
      corr=Matrix::bdiag(lapply(blocks,function(B){B[abs(B)<0.02]=0; methods::as(B,"dgCMatrix")}))
      sf=bigsparser::as_SFBM(methods::as(corr,"dgCMatrix"), backingfile=tempfile())
      bt=tryCatch({
            a=bigsnpr::snp_ldpred2_auto(sf,df,h2_init=0.3,vec_p_init=c(0.01,0.1),
                                        burn_in=100,num_iter=100,ncores=1)
            bb=rowMeans(sapply(a,function(z) z$beta_est)); if(all(is.na(bb))) stop("auto diverged"); bb
          }, error=function(e) bigsnpr::snp_ldpred2_inf(sf,df,h2=0.3))
      as.numeric(Ate %*% bt)
    })
    add(k,"LDpred2",r,if(r$ok) r$v else NULL,yts)
  }

  # ---- glinternet ---------------------------------------------------------
  if (has("glinternet") && want("glinternet") && p <= P_MAX_GLINTERNET) {
    r = budgeted({
      # screenLimit caps how many main effects can spawn interaction candidates;
      # without it the p^2 search dominates and never returns.
      fit=glinternet::glinternet.cv(Atr,ytr,numLevels=rep(1,p),nFolds=NFOLDS,
                                    nLambda=20,screenLimit=SCREENLIMIT)
      as.numeric(predict(fit,Ate,type="response"))
    })
    add(k,"glinternet",r,if(r$ok) r$v else NULL,yts)
  } else if (has("glinternet") && want("glinternet")) {
    add(k,"glinternet",list(v=NULL,s=NA_real_,ok=FALSE),NULL,yts,
        sprintf("infeasible: p=%d exceeds %d", p, P_MAX_GLINTERNET))
  }

  # ---- hierNet (SHIM stand-in) --------------------------------------------
  if (has("hierNet") && want("hiernet") && p <= P_MAX_HIERNET) {
    r = budgeted({
      # hierNet.cv multiplies cost by nfolds and was the dominant term; take a
      # fixed mid-path lambda instead so the fit completes.
      fit=hierNet::hierNet.path(Atr,ytr,nlam=5,trace=0)
      as.numeric(stats::predict(fit, newx=Ate)[,3])
    })
    add(k,"hierNet",r,if(r$ok) r$v else NULL,yts)
  } else if (has("hierNet") && want("hiernet")) {
    add(k,"hierNet",list(v=NULL,s=NA_real_,ok=FALSE),NULL,yts,
        sprintf("infeasible: p=%d exceeds %d", p, P_MAX_HIERNET))
  }

  # ---- shim, only if the 2017 GitHub build installed ----------------------
  if (has("shim") && want("shim")) {
    r = budgeted({ fit=shim::shim(Atr,ytr); as.numeric(Ate %*% as.numeric(stats::coef(fit))) })
    add(k,"shim",r,if(r$ok) r$v else NULL,yts)
  }

  rm(Atr,Ate,Xtr,Xte); gc(full=TRUE)
  if(k%%5==0){cat(sprintf("  rep %3d/%d\n",k,REPS));flush.console()}
}
res=do.call(rbind,rows)
out=paste0(RT,"result/sim",TAG,"_extra"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
OUTTAG = Sys.getenv("OUTTAG","run")
write.csv(res,sprintf("%s/eval_%s_%d_%s_%s.csv",out,simfun,nSNP,format(causal),OUTTAG),row.names=FALSE)
cat("[",id,"] DONE rows=",nrow(res),"\n",sep="")
