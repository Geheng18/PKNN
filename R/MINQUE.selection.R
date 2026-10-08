#' R-Language MINQUE
#'
#' Given a model specified by a response variable \code{y}, a series
#' of K kernels \code{V_1 ... V_K}, and an optional matrix of
#' covariate, solve the variance components and fixed coefficients.
#'
#' @param y. vector of response variable
#' @param v. list of kernel matrices
#' @param x. matrix of covariate
#' @param w. initial values for variance components
#' @return list of results, includs:
#' * vcs: vector of variance component estimates
#' * fix: vector of fixed coefficients, if the matrix covariate
#' \code{x.} was not null.
# ', ^ indicate transpose and generalized inverse; 
# R V_i and R y indicate Rv[[i]], Ry respectively;
# W is an initial value of VCs in redidual VC(1), genetic VC(2), ..., gecetic VC(K) order;
# rln_mnq corresponds to .mnq() in mnq.R file, with a little altered;
# 1 Difference between rln_mnq/_ginv() and .mnq():
#   V^ is by chol2inv(chol(V)) in rln_mnq() or by ginv(V) in rln_mnq(),
#   V^ is by solve(V, *) in .mnq()
#   the above difference appears when calcualting Rv[[i]], Ry
# 2 Some tricks are both used in rln_mnq/_ginv() and .mnq(), i.e., tric 1 and trick 2
# 3 The variance component algorithms used in LMM_MINQUE_Solver() in minque_solver.R file and in rln_mnq/_ginv(),
#   .mnq() are different, the LMM_MINQUE_Solver() follows the paper of Xiaoxi
#   type= 'MINQUE0' 'MINQUE1'
MINQUE.selection = function(y,geno,input.kernel,MINQUE.type = 'MINQUE1',KNN.type = 'interaction',lambda = 500,w = NULL,X = NULL,constrain = FALSE) {
  
  v = KNN2LMM(geno, input.kernel, KNN.type)
  K = length(v)
  N = length(y)
  
  
  if (MINQUE.type == 'MINQUE0') {
    A = diag(1, N, N)
  }
  
  if (MINQUE.type == 'MINQUE1') {
    ## initial weights
    if (is.null(w))
      w = rep(1 / K, K)
    
    ## sum of V_i, i = 1 .. K by initial weights
    # mapply(`*`, v., w., SIMPLIFY=FALSE) returns a list with each element as kernel * (corresponding)initial weight;
    # V is a matrix, which is v.[[1]] * w.[1] + ... + v.[[K]] * w.[K];
    V = Reduce(`+`, mapply(`*`, v, w, SIMPLIFY = FALSE)) # V is symmetric
    # generalized inverse of V, which is V^; if V is not positive definite, chol() will fail;
    # same as LMM_MINQUE_Solver() in minque_solver.R;
    A = solve(V)
  }
  
  ## get P, Q = I - P, and R V_i and R y, where R = V^Q
  Rv = list()
  if (is.null(X))
  {
    ## no X, P = 0, Q = I, R V_i = V^ I V_i = V^V_i
    for (i in seq.int(K))
      Rv[[i]] = A %*% v[[i]]
    Ry = A %*% y
  }
  else
  {
    ## P = X (X'V^X)^ (V^X)'        # ', ^ indicate transpose and generalized inverse;
    B = solve(V, X)                # V^X, solve V %*% b = X, then ^ indicates generalized inverse;
    C = MASS::ginv(t(X) %*% B) %*% t(B) # (X'V^X)^ (V^X)'
    P = X %*% C
    
    ## R V_i = V^ (I - P) V_i = V^ (V_i - P V_i)
    for (i in seq.int(K))
      Rv[[i]] = A %*% (v[[i]] - P %*% v[[i]])
    Ry = A %*% (y - P %*% y)
  }
  
  ## u_i = e' V_i e = y' R V_i R y
  u = double(K)                        # set u as double presion vector;
  for (i in seq.int(K))
    u[i] = sum(Ry * v[[i]] %*% Ry) # trick 1
  
  ## Caculate F: F_ij = Tr(R V_i R V_j)
  F = matrix(.0, K, K)
  for (i in seq.int(K))
  {
    for (j in seq.int(K)[-1L])
    {
      F[i, j] = sum(Rv[[i]] * Rv[[j]]) # trick 2, tr(A%*%B)=sum(A*B), if A,B are symmetric
      F[j, i] = F[i, j]
    }
    F[i, i] = sum(Rv[[i]] * Rv[[i]])
  }
  penalty = diag(lambda, K, K)
  penalty[1, 1] = 0
  F = F + penalty
  Cmat = MASS::ginv(F)
  
  ## [s2_1, .., s2_K]' = [yA1y, .., yAKy]' = solve(F, u) = F^u
  w = Cmat %*% u
  ##  bypass A1, .. A_K in (Rao. 1971)
  
  ## GLS for fixed effects
  if (!is.null(X))
  {
    b = C %*% y
    names(b) = colnames(X)
  }
  else
    b = NULL
  
  ## pack up
  theta = as.numeric(w)
  beta = as.numeric(b)
  if (constrain) {
    theta[theta < 0] = 0
  }
  
  return(list(theta = theta, beta = beta, Cmat = Cmat, v=v))
}
