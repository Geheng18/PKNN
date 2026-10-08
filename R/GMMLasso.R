
#'  A Lasso-based linear mixed model with generalized method-of-moments estimation for complex phenotype prediction.
#'
#' @param Y The continuous phenotypes of subjects.
#' @param Gen Gen is a list where each element is a n*p genomic matrix with n subjects and p SNPs.
#' @param index index is the index for testing subjects.
#' @param K K is a list where each element is the kernel matrix of the corresponding genomic matrix in the list Gen.
#' @param returnK returnK=T means the kernel matrix will be returned.returnK=F means the kernel matrix will not be returned. The default is returnK=F
#'
#'
#' @return The prediction values and ture values for testing subjects will be returned.
#' @return The effects sizes for each genomic region can also returned.
#' @export
#' @import glmnet
#' @import MASS
#' @importFrom glmnet cv.glmnet
#'
#' @examples
#'OGen=matrix(sample(0:2,500*150,replace = TRUE),500,150)
#'start <- seq(1, by = 3, length = ncol(OGen) / 3)
#'Gen <- lapply(start, function(i, OGen) OGen[,i:(i+2)], OGen = OGen)
#'y = rowSums(scale(Gen[[1]]))+rowSums(scale(Gen[[2]]))+rnorm(500)
#'index= sample(1:length(y),100)
#'fit=GmmLasso::GMMLasso(y=y, Gen=Gen, index=index, K = NULL, returnK = F)

GMMLasso<-function(y,X=NULL,Gen,index,K=NULL,returnK=F,lambda=NULL){

    test_index <- index  ##test index
    train_index <- c(1:length(y))[-index] ##train index


    Pheno <- y[train_index]  ##training Phenotypes
    mean_adjust=mean(Pheno)
    sd_adjust=sd(Pheno)
    Pheno=scale(Pheno)



    if(is.null(K)){
      K=lapply(Gen,Fun_Klinear) ##get kernel function
    }
    list_K=lapply(K,function(kk){
      kk[train_index,train_index]
    })

    list_K=lapply(list_K, scaleK)


    numK=length(list_K)
    list_K[[numK+1]]=diag(length(train_index))

    if(is.null(X)){
      n=length(Pheno);
      V=diag(n)
      E=eigen(V)$vectors
      A=E                       ##get A matrix
    }else{
      n=length(Pheno);
      V=diag(n)-X %*% MASS::ginv(t(X) %*% X) %*% t(X)
      E=eigen(V)$vectors
      nc=ncol(X)
      A=E[,1:(n-nc)]                       ##get A matrix
    }


    lasso_Pheno = c(t(A)%*%Pheno%*%t(Pheno)%*%A)   ##training Phenotypes used in GMMLasso
    Tlist=lapply(list_K,function(kk) c(t(A)%*%kk%*%A))
    T=matrix(unlist(Tlist),nrow=ncol(A)*ncol(A),ncol=length(list_K))


    p.fac = c(rep(1,numK),0)
    fit_cv <- cv.glmnet(T, lasso_Pheno, alpha=0,lower.limits =0,penalty.factor = p.fac)     ## fit the model
	
	if(is.null(lambda)){
		lambda=chooselambda(fit_cv,1) 
	}
        ##search the optimal lambda
    

    fit_best <- glmnet(T, lasso_Pheno, alpha = 1, lambda = lambda, lower.limits =0, penalty.factor = p.fac)       ## fit the model
    beta= t(as.matrix(fit_best$beta))   #get the parameter estimations
    mu=coef(fit_best)[1]



    Sigma=matrix(0,length(y),length(y))

    for (j in 1:numK) {
      Sigma=Sigma + fit_best$beta[j] * K[[j]]
    }

    Sigma=Sigma+fit_best$beta[length(list_K)]*diag(length(y))           ##get the variance matrix


    Sigma11=Sigma[test_index,test_index]
    Sigma12=Sigma[test_index,train_index]
    Sigma21=Sigma[train_index,test_index]
    Sigma22=Sigma[train_index,train_index]



    Pheno_pred=mu+Sigma12 %*% MASS::ginv(Sigma22) %*%(Pheno-mu)
    Pheno_pred=Pheno_pred*sd_adjust+mean_adjust

    combineout=cbind(pred=Pheno_pred,true=y[test_index])

    if(returnK)  out=list(out=combineout,beta=beta, K=K)
    else out=list(out=combineout,beta=beta)
    out
  }

