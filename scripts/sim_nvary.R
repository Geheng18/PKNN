# ---------------------------------------------------------------------------
# Accuracy as a function of training sample size (Reviewer 2, comment 2).
#
# Fixed:  20 SNPs per gene, 10 regions, causal rate 0.2, h2 = 0.5.
#         All NINE architectures (5 link functions + 4 interaction mechanisms),
#         so the sweep shows whether the advantage of the proposed method grows
#         with sample size rather than only how a linear model converges.
# Varied: N = 500, 1000, 2000, 3000, 5000, 10000 -- the same grid as the runtime
#         scalability study, with the same split rule nTest = max(100, 0.2N), so
#         N = 1000 reproduces the 800/200 split of the main simulation exactly
#         and the two studies are directly comparable.
#
# Data are generated inside the loop rather than read from disk: at N = 10,000
# the genotypes for 100 replicates would be far larger than the results.
#
# TASK TABLE.  Cost per replicate was measured in the runtime study, and tasks
# are sized from it so none exceeds ~15 h and none exceeds its memory request.
# Large N is split by replicate range, and at N = 10,000 also by method, because
# the two numerical MINQUE variants need ~450 GB there while everything else
# fits in 360 GB.  Replicate counts fall with N for the same reason; the count
# actually used is recorded in every output row.
#
# Output: result/sim_nvary/eval_N<N>_<tag>_<from>-<to>.csv
# ---------------------------------------------------------------------------
suppressMessages({library(BEDMatrix); library(data.table); library(MASS)
                  library(glmnet); library(matrixcalc); library(Matrix)})
LIB="/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/Rlib"; .libPaths(c(LIB,.libPaths()))
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
for (f in list.files(paste0(RT,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
UKB = "/home/heng.ge/blue/heng.ge/KNN_Data_backup_again/data/UKB/bchQC03"
has = function(p) requireNamespace(p, quietly=TRUE)

nSNP = 20; nRegion = 10; nRegion.true = 2; H2 = 0.5
NONLIN = c("linear","cosh","square","hyperbola","rickercurve")
INTER  = c("within.region.threshold.TRUE","within.region.threshold.FALSE",
           "outside.region.threshold.TRUE","outside.region.threshold.FALSE")
SCENS  = c(NONLIN, INTER)
ALL  = c("CF-ridge","NUM-lasso","NUM-enet","GMM","BLUP","LDpred2","hierNet")
CHEAP = c("CF-ridge","GMM","BLUP","LDpred2","hierNet")   # all but the NUM pair
NUMS  = c("NUM-lasso","NUM-enet")

# Chunks per architecture, sized from the per-replicate costs measured in the
# runtime study so no task exceeds ~14 h.  N = 10,000 is omitted here: it needs
# 400-540 GB, and by N = 5,000 the methods have largely converged, so the extra
# level buys little for nine architectures.  The linear architecture was already
# run over the full grid, including N = 10,000, in a first pass.
CHUNKS = rbind(
  data.frame(N=  500, from= 1, to=100),
  data.frame(N= 1000, from= 1, to=100),
  data.frame(N= 2000, from= 1, to=100),
  data.frame(N= 3000, from= 1, to= 50),
  data.frame(N= 3000, from=51, to=100),
  data.frame(N= 5000, from= 1, to= 20),
  data.frame(N= 5000, from=21, to= 40),
  data.frame(N= 5000, from=41, to= 60),
  data.frame(N= 5000, from=61, to= 80),
  data.frame(N= 5000, from=81, to=100))
# scenario varies slowest: indices (s-1)*10 + 1 .. s*10 belong to architecture s
TASKS = do.call(rbind, lapply(SCENS, function(sc)
          cbind(CHUNKS, scen=sc, tag="all", methods=paste(ALL,collapse="|"),
                stringsAsFactors=FALSE)))

id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); stopifnot(id>=1, id<=nrow(TASKS))
tk = TASKS[id,]; N=tk$N; scen=tk$scen; want=strsplit(tk$methods,"\\|")[[1]]
# NVARY_FROM / NVARY_TO override the chunk's replicate range, for topping up a
# chunk that hit the wall clock.  Output is named from the range, so a top-up
# writes its own file rather than overwriting the partial one.
if (nzchar(Sys.getenv("NVARY_FROM"))) tk$from = as.integer(Sys.getenv("NVARY_FROM"))
if (nzchar(Sys.getenv("NVARY_TO")))   tk$to   = as.integer(Sys.getenv("NVARY_TO"))
is.inter = scen %in% INTER
nTest = max(100, round(0.2*N)); nTrain = N - nTest
cat(sprintf("[task %d/%d] %s  N=%d nTrain=%d reps %d-%d  methods: %s\n",
            id, nrow(TASKS), scen, N, nTrain, tk$from, tk$to, paste(want,collapse=",")))
flush.console()

bm  = BEDMatrix(paste0(UKB,".bed"), simple_names=TRUE)
bim = fread(paste0(UKB,".bim"), header=FALSE)
cs  = function(x) scale(x, scale=FALSE)
tm  = function(e){t0=Sys.time()
  v=tryCatch(force(e), error=function(z) structure(list(m=conditionMessage(z)),class="fx"))
  list(v=v, s=as.numeric(difftime(Sys.time(),t0,units="secs")), ok=!inherits(v,"fx"),
       msg=if(inherits(v,"fx")) v$m else NA_character_)}
rows=list()
add=function(rep,method,r,pred,yts){
  mse=NA_real_; cc=NA_real_
  if (r$ok && !is.null(pred) && all(is.finite(pred)) && sd(pred)>0){
    mse=mean((yts-pred)^2); cc=cor(yts,pred) }
  rows[[length(rows)+1]] <<- data.frame(rep=rep, scen=scen, N=N, nTrain=nTrain,
    method=method, mse=mse, cor=cc, secs=r$s, ok=r$ok, msg=r$msg) }
flush.rows=function(){
  if(!length(rows)) return(invisible())
  dir.create(paste0(RT,"result/sim_nvary"), recursive=TRUE, showWarnings=FALSE)
  write.csv(do.call(rbind,rows), sprintf("%sresult/sim_nvary/eval_%s_N%d_%s_%d-%d.csv",
            RT, scen, N, tk$tag, tk$from, tk$to), row.names=FALSE) }

for (k in tk$from:tk$to) {
  set.seed(20261005 + 10000*N + 7919*match(scen,SCENS) + k)                 # same data for every method
  gt = tm(if (is.inter)
            Simdata.LD.INT(bm, bim, N=N, nSNP=nSNP, nRegion=nRegion,
                           nRegion.true=nRegion.true, scenario=scen, h2=H2)
          else
            Simdata.LD(bm, bim, N=N, nSNP=nSNP, nRegion=nRegion,
                       nRegion.true=nRegion.true, simfun=scen, h2=H2))
  if (!gt$ok) { cat("  rep",k,"generation failed:",gt$msg,"\n"); next }
  d = gt$v
  Xtr  = lapply(d$X, function(G) cs(G[1:nTrain,]))
  Xall = lapply(d$X, function(G) rbind(cs(G[1:nTrain,]), cs(G[(nTrain+1):N,])))
  ytr = scale(d$y[1:nTrain]); yts = scale(d$y[(nTrain+1):N])
  ik  = rep(list(c("product")), nRegion); idx=(nTrain+1):N
  cat(sprintf("  rep %d  (generated in %.0fs)\n", k, gt$s)); flush.console()

  if ("CF-ridge" %in% want) {
    r = tm(MINQUE.selection(ytr, Xtr, ik, MINQUE.type="MINQUE1", lambda=500, constrain=TRUE))
    p = if(r$ok) tryCatch(MINQUE.predict(ytr,Xall,ik,r$v$theta), error=function(e) NULL) else NULL
    add(k,"CLOSEFORM-MINQUE1-RIDGE",r,p,yts); flush.rows() }
  for (a in list(c("NUM-lasso","NUMERICAL-MINQUE-LASSO",1),
                 c("NUM-enet","NUMERICAL-MINQUE-ELASTICNET",0.5))) {
    if (!(a[1] %in% want)) next
    r = tm(MINQUE.selection.numerical(ytr, Xtr, ik, alpha=as.numeric(a[3])))
    p = if(r$ok) tryCatch(MINQUE.predict(ytr,Xall,ik,r$v$theta), error=function(e) NULL) else NULL
    add(k,a[2],r,p,yts); flush.rows() }

  if (any(c("GMM","BLUP") %in% want)) {
    K = lapply(Xall, function(G) findKernel("product", G)); yf=c(ytr,yts)
    if ("GMM" %in% want) {
      r = tm(GMMLasso(y=yf, index=idx, K=K))
      add(k,"GMM",r,if(r$ok) as.numeric(r$v$out[,1]) else NULL,yts); flush.rows() }
    if ("BLUP" %in% want && has("rrBLUP")) {
      r = tm({ Ka=scaleK(Reduce(`+`,K)); Ktt=Ka[1:nTrain,1:nTrain]; Kts=Ka[idx,1:nTrain]
               ms=rrBLUP::mixed.solve(y=as.numeric(ytr), K=Ktt)
               as.numeric(Kts %*% solve(Ktt+(ms$Ve/ms$Vu)*diag(nTrain), as.numeric(ytr))) })
      add(k,"BLUP",r,if(r$ok) r$v else NULL,yts); flush.rows() }
    rm(K); gc(full=TRUE) }

  if (any(c("LDpred2","hierNet") %in% want)) {
    Atr=do.call(cbind,Xtr); Ate=do.call(cbind,lapply(d$X,function(G) cs(G[(nTrain+1):N,])))
    sdv=apply(Atr,2,sd); keep=which(sdv>0); Atr=Atr[,keep,drop=FALSE]; Ate=Ate[,keep,drop=FALSE]
    if ("LDpred2" %in% want && has("bigsnpr")) {
      r = tm({
        b=as.numeric(crossprod(Atr,ytr)/colSums(Atr^2))
        res=as.numeric(ytr) - sweep(Atr,2,b,`*`)
        se=sqrt(colSums(res^2)/(nTrain-2)/colSums(Atr^2))
        df=data.frame(beta=b, beta_se=se, n_eff=nTrain)
        off=0; blocks=list()
        for (l in 1:nRegion) { m=ncol(Xtr[[l]]); jl=intersect(keep,(off+1):(off+m))-off
          if(length(jl)>1) blocks[[length(blocks)+1]]=cor(Xtr[[l]][,jl,drop=FALSE]); off=off+m }
        corr=Matrix::bdiag(lapply(blocks,function(B){B[abs(B)<0.02]=0; methods::as(B,"dgCMatrix")}))
        sf=bigsparser::as_SFBM(methods::as(corr,"dgCMatrix"), backingfile=tempfile())
        bt=tryCatch({ au=bigsnpr::snp_ldpred2_auto(sf,df,h2_init=0.3,vec_p_init=c(0.01,0.1),
                        burn_in=100,num_iter=100,ncores=1)
                      bb=rowMeans(sapply(au,function(z) z$beta_est))
                      if(all(is.na(bb))) stop("auto diverged"); bb },
                    error=function(e) bigsnpr::snp_ldpred2_inf(sf,df,h2=0.3))
        as.numeric(Ate %*% bt) })
      add(k,"LDpred2",r,if(r$ok) r$v else NULL,yts); flush.rows() }
    if ("hierNet" %in% want && has("hierNet")) {
      r = tm({ fit=hierNet::hierNet.path(Atr, as.numeric(ytr), nlam=5, trace=0)
               as.numeric(stats::predict(fit, newx=Ate)[,3]) })
      add(k,"hierNet",r,if(r$ok) r$v else NULL,yts); flush.rows() }
    rm(Atr,Ate); gc(full=TRUE) }

  rm(d,Xtr,Xall); gc(full=TRUE)
}
flush.rows()
cat("[task ",id,"] DONE rows=",length(rows),"\n",sep="")
