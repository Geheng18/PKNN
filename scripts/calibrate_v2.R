# Calibration for the corrected design: fixed number of causal SNPs per causal
# region, so increasing nSNP adds noise SNPs instead of diluting the signal.
suppressMessages({library(BEDMatrix); library(MASS); library(glmnet); library(data.table)})
RT="/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
UKB="/home/heng.ge/blue/heng.ge/KNN_Data_backup_again/data/UKB/bchQC03"
for (f in list.files(paste0(RT,"R"),pattern="\\.R$",full.names=TRUE)) source(f)
grid = expand.grid(simfun=c("linear","cosh","rickercurve"), nSNP=c(20,100,500,1000),
                   h2=c(0.3,0.5), ncs=20, stringsAsFactors=FALSE)
id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID","1")); cf=grid[id,]
N=1000; nTrain=800; nRegion=10; REPS=5
bm=BEDMatrix(paste0(UKB,".bed"),simple_names=TRUE); bim=fread(paste0(UKB,".bim"),header=FALSE)
cs=c()
for (rep in 1:REPS) {
  set.seed(id*991L+rep)
  d=Simdata.LD(bm,bim,N=N,nSNP=cf$nSNP,nRegion=nRegion,nRegion.true=4,
               simfun=cf$simfun,h2=cf$h2,n.causal.snp=cf$ncs)
  Xtr=lapply(d$X,function(G) scale(G[1:nTrain,],scale=FALSE))
  Xall=lapply(d$X,function(G) rbind(scale(G[1:nTrain,],scale=FALSE),scale(G[(nTrain+1):N,],scale=FALSE)))
  ytr=scale(d$y[1:nTrain]); yts=scale(d$y[(nTrain+1):N]); ik=rep(list(c('product')),nRegion)
  fit=MINQUE.selection(ytr,Xtr,ik,MINQUE.type='MINQUE1',lambda=500,constrain=TRUE)
  cs=c(cs, cor(yts,MINQUE.predict(ytr,Xall,ik,fit$theta)))
  cat("rep",rep,"cor",round(tail(cs,1),3),"\n"); flush.console()
}
res=data.frame(simfun=cf$simfun,nSNP=cf$nSNP,h2=cf$h2,n.causal.snp=cf$ncs,
               meanCOR=mean(cs,na.rm=TRUE),sdCOR=sd(cs,na.rm=TRUE))
dir.create(paste0(RT,"result/calib2"),recursive=TRUE,showWarnings=FALSE)
write.csv(res,sprintf("%sresult/calib2/c2_%s_P%d_h%.2f.csv",RT,cf$simfun,cf$nSNP,cf$h2),row.names=FALSE)
print(res)
