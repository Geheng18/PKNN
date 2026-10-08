# Timing probe: how expensive is one replicate at each region count?
# The interaction basis has 2 + R + R(R+1)/2 components, so cost grows fast.
suppressMessages({library(MASS); library(glmnet); library(matrixcalc)})
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
for (f in list.files(paste0(RT,"R"), pattern="\\.R$", full.names=TRUE)) source(f)
nTrain=800; nTest=200; N=nTrain+nTest
center_scale = function(x) scale(x, scale=FALSE)
tm = function(e){t0=Sys.time(); v=tryCatch(force(e),error=function(z)NULL)
                 list(v=v,s=as.numeric(difftime(Sys.time(),t0,units="secs")))}
for (R in as.integer(strsplit(Sys.getenv("PROBE_R","10+20"),"[,+]")[[1]])) {
  cell = sprintf("%sdata/sim_ngene/linear/%d/", RT, R)
  if (!file.exists(sprintf("%sSim_1_1_%d.Rdata", cell, R))) { cat("R=",R," no data\n"); next }
  Xtr=list(); Xall=list()
  for (l in 1:R) { load(sprintf("%sSim_1_%d_%d.Rdata", cell, l, R)); G=as.matrix(geno)
    Xtr[[l]]=apply(G[1:nTrain,],2,center_scale)
    Xall[[l]]=rbind(apply(G[1:nTrain,],2,center_scale),apply(G[(nTrain+1):N,],2,center_scale)) }
  y=unlist(read.table(sprintf("%sPhe_1.txt",cell))); ytr=scale(y[1:nTrain]); yts=scale(y[(nTrain+1):N])
  ik=rep(list(c("product")),R)
  nComp = 2 + R + R*(R+1)/2
  a=tm(MINQUE.selection(ytr,Xtr,ik,MINQUE.type="MINQUE1",lambda=500,constrain=TRUE))
  b=tm(if(!is.null(a$v)) MINQUE.predict(ytr,Xall,ik,a$v$theta) else NULL)
  d=tm(MINQUE.selection.numerical(ytr,Xtr,ik,alpha=1))
  K=lapply(Xall,function(G) findKernel("product",G))
  g=tm(GMMLasso(y=c(ytr,yts), index=(nTrain+1):N, K=K))
  cat(sprintf("R=%-3d nComp=%-4.0f CF=%7.1fs pred=%6.1fs NUMlasso=%8.1fs GMM=%7.1fs  peakGB=%.1f\n",
      R, nComp, a$s, b$s, d$s, g$s, sum(gc()[,6])/1024)); flush.console()
}
