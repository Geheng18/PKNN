# ---------------------------------------------------------------------------
# "How many candidate genes?" experiment -- reviewer-requested comparators.
#   LDpred2  (bigsnpr)  and  hierNet  (explicit pairwise interactions)
#
# Same 40-cell grid as analyze_ngene.R.  Both comparators are feasible at every
# region count here: p = nRegion * 20 reaches only 400 at R = 20, below the
# hierNet ceiling of 500 established in the main simulation.
#
# Output: result/sim_ngene_extra/eval_<simfun>_<R>_<method>.csv
# ---------------------------------------------------------------------------
LIB="/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/Rlib"; .libPaths(c(LIB,.libPaths()))
suppressMessages({library(MASS); library(Matrix)})
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
has = function(p) requireNamespace(p, quietly=TRUE)
BUDGET = as.numeric(Sys.getenv("BUDGET","900"))
P_MAX_HIERNET = as.numeric(Sys.getenv("P_MAX_HIERNET","500"))
WANT = strsplit(Sys.getenv("METHODS","ldpred2+hiernet"), "[,+]")[[1]]
want = function(m) m %in% WANT

nTrain=800; nTest=200; N=nTrain+nTest
center_scale = function(x) scale(x, scale=FALSE)

SCEN   = strsplit(Sys.getenv("SCENARIOS","linear+square+cosh+hyperbola+rickercurve"), "[,+]")[[1]]
REGION = as.integer(strsplit(Sys.getenv("REGIONS","2+4+6+8+10+12+16+20"), "[,+]")[[1]])
grid = expand.grid(simfun=SCEN, nRegion=REGION, stringsAsFactors=FALSE)
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id>=1, id<=nrow(grid))
cf = grid[id,]; simfun=cf$simfun; nRegion=cf$nRegion
REPS = as.integer(Sys.getenv("REPS","100"))
cell = sprintf("%sdata/sim_ngene/%s/%d/", RT, simfun, nRegion)
cat(sprintf("[%d] %s nRegion=%d p=%d reps=%d methods=%s\n", id, simfun, nRegion,
            nRegion*20, REPS, paste(WANT,collapse=","))); flush.console()

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
flush.rows=function(){
  if(!length(rows)) return(invisible())
  d=do.call(rbind,rows); dir.create(paste0(RT,"result/sim_ngene_extra"),recursive=TRUE,showWarnings=FALSE)
  for (m in unique(d$method)) {
    tag = if (m=="LDpred2") "ldpred2" else "hiernet"
    write.csv(d[d$method==m,], sprintf("%sresult/sim_ngene_extra/eval_%s_%d_%s.csv",
              RT, simfun, nRegion, tag), row.names=FALSE)
  }
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
  Atr=do.call(cbind,Xtr); Ate=do.call(cbind,Xte)
  y=unlist(read.table(sprintf("%sPhe_%d.txt",cell,k)))
  ytr=as.numeric(scale(y[1:nTrain])); yts=as.numeric(scale(y[(nTrain+1):N]))
  sdv=apply(Atr,2,stats::sd); keep=which(sdv>0)
  Atr=Atr[,keep,drop=FALSE]; Ate=Ate[,keep,drop=FALSE]; p=ncol(Atr)

  if (has("bigsnpr") && want("ldpred2")) {
    r = budgeted({
      n=nTrain
      b=as.numeric(crossprod(Atr,ytr)/colSums(Atr^2))
      res=ytr - sweep(Atr,2,b,`*`)
      se=sqrt(colSums(res^2)/(n-2)/colSums(Atr^2))
      df=data.frame(beta=b, beta_se=se, n_eff=n)
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

  if (has("hierNet") && want("hiernet") && p <= P_MAX_HIERNET) {
    r = budgeted({
      fit=hierNet::hierNet.path(Atr,ytr,nlam=5,trace=0)
      as.numeric(stats::predict(fit, newx=Ate)[,3])
    })
    add(k,"hierNet",r,if(r$ok) r$v else NULL,yts)
  } else if (has("hierNet") && want("hiernet")) {
    add(k,"hierNet",list(v=NULL,s=NA_real_,ok=FALSE),NULL,yts,
        sprintf("infeasible: p=%d exceeds %d", p, P_MAX_HIERNET))
  }

  rm(Atr,Ate,Xtr,Xte); gc(full=TRUE)
  if (k %% 5 == 0) { flush.rows(); cat(sprintf("  rep %3d/%d\n", k, REPS)); flush.console() }
}
flush.rows()
cat("[",id,"] DONE rows=",length(rows),"\n",sep="")
