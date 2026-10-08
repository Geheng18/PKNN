# ---------------------------------------------------------------------------
# Task 2 -- runtime and memory vs sample size (N) and SNP count (P)
#
# Answers Reviewer 1 (Model Design Q2): "In addition to runtime, I would also
# suggest reporting the memory usage with increasing sample size, as well as
# number of SNPs".
#
# Two bases are benchmarked, because the memory wall is set by the variance-
# component basis, not by MINQUE itself:
#   interaction : I, J, {J*K_i}, {K_i*K_j}  -> 67 components for 10 regions
#   first       : I, J, K_1..K_10           -> 12 components
#
# Each array task = one (basis, N, nSNP) config, run as its own SLURM job so
# job-level MaxRSS from sacct is a clean per-config peak-memory measurement.
# Methods that exceed memory/time are recorded as failures -- that boundary is
# itself the scalability result.
# ---------------------------------------------------------------------------
suppressMessages({library(BEDMatrix); library(MASS); library(glmnet); library(data.table)})

RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
UKB = "/home/heng.ge/blue/heng.ge/KNN_Data_backup_again/data/UKB/bchQC03"
for (f in list.files(paste0(RT, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

nRegion = 10; H2 = 0.20; REPS = 3

# All benchmark cells use the interaction basis (67 components for 10 regions),
# which is the basis the method actually uses.  The earlier first-order rows were
# retired; the N=3000 ceiling there was an allocation limit of mine (160 GB), not
# a property of the basis -- 67 dense n x n components need ~100 GB at N=10,000.
grid = rbind(
  expand.grid(basis = "interaction", N = c(500,1000,2000,3000), nSNP = c(20,100,500,1000), stringsAsFactors = FALSE),  #  1-16 done
  expand.grid(basis = "interaction", N = c(5000,10000),         nSNP = c(20,100,500,1000), stringsAsFactors = FALSE),  # 17-24
  expand.grid(basis = "interaction", N = c(5000,10000),         nSNP = c(5000,10000),      stringsAsFactors = FALSE),  # 25-28
  expand.grid(basis = "interaction", N = c(500,1000,2000,3000), nSNP = c(5000,10000),      stringsAsFactors = FALSE))  # 29-36

id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot(id >= 1, id <= nrow(grid))
cf = grid[id, ]
basis = cf$basis; N = cf$N; nSNP = cf$nSNP
nTest = max(100, round(0.2 * N)); nTrain = N - nTest

cat(sprintf("[task %d] basis=%s N=%d nSNP=%d nTrain=%d\n", id, basis, N, nSNP, nTrain))
flush.console()

bm  = BEDMatrix(paste0(UKB, ".bed"), simple_names = TRUE)
bim = fread(paste0(UKB, ".bim"), header = FALSE)

# peak R heap (MB) around an expression, plus wall time and any error
timed = function(expr) {
  gc(reset = TRUE, full = TRUE)
  t0 = Sys.time()
  val = tryCatch(force(expr), error = function(e) structure(list(msg = conditionMessage(e)), class = "bmfail"))
  el  = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  g   = gc(full = TRUE)
  list(value = val, secs = el, peakMB = sum(g[, 6]),
       ok = !inherits(val, "bmfail"),
       msg = if (inherits(val, "bmfail")) val$msg else NA_character_)
}

WANT = strsplit(Sys.getenv("BENCH_METHODS","cf+gmm+numlasso+numenet+blup+ldpred2+hiernet"), "[,+]")[[1]]
want = function(x) x %in% WANT
SUF  = Sys.getenv("BENCH_SUFFIX","")
HN.BUDGET = as.numeric(Sys.getenv("HIERNET_BUDGET","10800"))   # seconds per hierNet fit
.libPaths(c(file.path(RT,"Rlib"), .libPaths()))
has = function(q) requireNamespace(q, quietly=TRUE)
OUTFILE = file.path(RT, "result", "benchmark", sprintf("bench_%s_N%d_P%d%s.csv", basis, N, nSNP, SUF))
dir.create(dirname(OUTFILE), recursive = TRUE, showWarnings = FALSE)
# Written after every method, so an out-of-memory kill in one method no longer
# discards the methods that already finished in the same cell.
flush.rows = function() write.csv(do.call(rbind, rows), OUTFILE, row.names = FALSE)
MEMCAP = as.numeric(Sys.getenv("MEMCAP_GB", "160"))
nComp  = if (basis == "interaction") 1 + (nRegion+1)*(nRegion+2)/2 else nRegion + 2
# Peak-memory model for NUM-MINQUE-LASSO, fitted to measured cells
#   GB ~ 439 * (8 n^2 / 1e9) + 7.53 * (8 n^2 K / 1e9),  n = training size, K = components
# reproduces the N=2000/3000 interaction and N=5000/10000 first-order peaks to ~10%.
# It predicts ~480 GB at N=10000 on the 67-component basis.
num.est.GB = (439 * 8 * nTrain^2 + 7.53 * 8 * nTrain^2 * nComp) / 1e9
rows = list()
for (rep in seq_len(REPS)) {
  set.seed(id * 1000L + rep)
  d = Simdata.LD(bm, bim, N = N, nSNP = nSNP, nRegion = nRegion,
                 nRegion.true = 4, simfun = "cosh", h2 = H2)
  Xtr  = lapply(d$X, function(G) scale(G[1:nTrain, ], scale = FALSE))
  Xall = lapply(d$X, function(G) rbind(scale(G[1:nTrain, ], scale = FALSE),
                                       scale(G[(nTrain+1):N, ], scale = FALSE)))
  ytr = scale(d$y[1:nTrain]); yts = scale(d$y[(nTrain+1):N])
  ik  = rep(list(c("product")), nRegion)

  # --- closed-form MINQUE, ridge penalty -----------------------------------
  if (want("cf")) {
  r1 = timed(MINQUE.selection(ytr, Xtr, ik, MINQUE.type = "MINQUE1", KNN.type = basis,
                              lambda = 500, constrain = TRUE))
  c1 = if (r1$ok) tryCatch(cor(yts, MINQUE.predict(ytr, Xall, ik, r1$value$theta, KNN.type = basis)),
                           error = function(e) NA_real_) else NA_real_
  rows[[length(rows)+1]] = data.frame(basis, N, nSNP, rep, method = "CF-MINQUE1-RIDGE",
                                      secs = r1$secs, peakMB = r1$peakMB, cor = c1,
                                      ok = r1$ok, msg = r1$msg)
  }
  flush.rows()
  # --- GMMLasso (region kernels; independent of the MINQUE basis) -----------------
  if (want("gmm")) {
    K.list = lapply(Xall, function(G) findKernel("product", G))
    r3 = timed(GMMLasso(y = c(ytr, yts), index = (nTrain+1):N, K = K.list))
    c3 = if (r3$ok) tryCatch(cor(yts, as.numeric(r3$value$out[,1])), error = function(e) NA_real_) else NA_real_
    rows[[length(rows)+1]] = data.frame(basis, N, nSNP, rep, method = "GMM",
                                        secs = r3$secs, peakMB = r3$peakMB, cor = c3,
                                        ok = r3$ok, msg = r3$msg)
    rm(K.list)
  }
  flush.rows()
  # --- numerical MINQUE, lasso ---------------------------------------------
  if (want("numlasso")) {
  r2 = if (num.est.GB > 0.85 * MEMCAP)
         list(value = NULL, secs = NA_real_, peakMB = NA_real_, ok = FALSE,
              msg = sprintf("infeasible: estimated %.0f GB exceeds %.0f GB allocation", num.est.GB, MEMCAP))
       else timed(MINQUE.selection.numerical(ytr, Xtr, ik, alpha = 1, KNN.type = basis))
  c2 = if (r2$ok) tryCatch(cor(yts, MINQUE.predict(ytr, Xall, ik, r2$value$theta, KNN.type = basis)),
                           error = function(e) NA_real_) else NA_real_
  rows[[length(rows)+1]] = data.frame(basis, N, nSNP, rep, method = "NUM-MINQUE-LASSO",
                                      secs = r2$secs, peakMB = r2$peakMB, cor = c2,
                                      ok = r2$ok, msg = r2$msg)
  }

  P = nRegion * nSNP                      # explicit predictors
  add.row = function(meth, r, cc = NA_real_)
    rows[[length(rows)+1]] <<- data.frame(basis, N, nSNP, rep, method = meth,
      secs = r$secs, peakMB = r$peakMB, cor = cc, ok = r$ok, msg = r$msg)

  # --- numerical MINQUE, elastic net ---------------------------------------
  if (want("numenet")) {
    r = if (num.est.GB > 0.85 * MEMCAP)
          list(value=NULL, secs=NA_real_, peakMB=NA_real_, ok=FALSE,
               msg=sprintf("infeasible: estimated %.0f GB exceeds %.0f GB allocation", num.est.GB, MEMCAP))
        else timed(MINQUE.selection.numerical(ytr, Xtr, ik, alpha = 0.5, KNN.type = basis))
    cc = if (r$ok) tryCatch(cor(yts, MINQUE.predict(ytr, Xall, ik, r$value$theta, KNN.type = basis)),
                            error=function(e) NA_real_) else NA_real_
    add.row("NUM-MINQUE-ENET", r, cc); flush.rows()
  }
  # --- linear BLUP on the aggregate GRM ------------------------------------
  if (want("blup") && has("rrBLUP")) {
    r = timed({ Kb = scaleK(Reduce(`+`, lapply(Xall, function(G) findKernel("product", G))))
                ms = rrBLUP::mixed.solve(y = as.numeric(ytr), K = Kb[1:nTrain,1:nTrain])
                lam = ms$Ve/ms$Vu
                as.numeric(Kb[(nTrain+1):N,1:nTrain] %*%
                           solve(Kb[1:nTrain,1:nTrain] + lam*diag(nTrain), as.numeric(ytr))) })
    cc = if (r$ok) tryCatch(cor(yts, r$value), error=function(e) NA_real_) else NA_real_
    add.row("BLUP", r, cc); flush.rows()
  }
  # --- LDpred2 --------------------------------------------------------------
  if (want("ldpred2") && has("bigsnpr")) {
    r = timed({
      G = do.call(cbind, Xall); sdv = apply(G[1:nTrain,,drop=FALSE], 2, sd)
      okc = is.finite(sdv) & sdv > 0
      Xa = sweep(G[1:nTrain,,drop=FALSE], 2, ifelse(okc, sdv, 1), `/`)
      sxx = colSums(Xa^2); xty = as.numeric(crossprod(Xa, ytr))
      b = ifelse(okc, xty/sxx, 0)
      se = sqrt(pmax(sum(ytr^2) - b*xty, 0)/(nTrain-2)/sxx)
      bad = !okc | !is.finite(se) | se <= 0; b[bad] = 0; se[bad] = 1
      blk = split(seq_len(ncol(G)), rep(seq_len(nRegion), times = sapply(Xall, ncol)))
      ii = vector("list", length(blk)); jj = ii; xx = ii
      for (t in seq_along(blk)) { bb = blk[[t]]
        C = suppressWarnings(stats::cor(Xa[, bb, drop=FALSE])); C[is.na(C)] = 0
        C[abs(C) < 0.02] = 0; diag(C) = 1
        nz = which(C != 0, arr.ind = TRUE)
        ii[[t]] = bb[nz[,1]]; jj[[t]] = bb[nz[,2]]; xx[[t]] = C[nz] }
      sf = bigsparser::as_SFBM(methods::as(Matrix::sparseMatrix(i=unlist(ii), j=unlist(jj),
             x=unlist(xx), dims=c(ncol(G),ncol(G))), "dgCMatrix"), backingfile=tempfile())
      df = data.frame(beta=b, beta_se=se, n_eff=nTrain)
      bt = tryCatch({ au = bigsnpr::snp_ldpred2_auto(sf, df, h2_init=0.2, vec_p_init=c(0.01,0.1),
                                                     burn_in=100, num_iter=100, ncores=1)
                      Bm = sapply(au, function(z) z$beta_est)
                      kc = apply(Bm, 2, function(z) all(is.finite(z)))
                      if (!any(kc)) stop("all chains diverged"); rowMeans(Bm[, kc, drop=FALSE]) },
                    error = function(e) bigsnpr::snp_ldpred2_inf(sf, df, h2=0.2))
      bt[!is.finite(bt)] = 0
      Xb = sweep(G[(nTrain+1):N,,drop=FALSE], 2, ifelse(okc, sdv, 1), `/`); Xb[, !okc] = 0
      as.numeric(Xb %*% bt) })
    cc = if (r$ok) tryCatch(cor(yts, r$value), error=function(e) NA_real_) else NA_real_
    add.row("LDpred2", r, cc); flush.rows()
  }
  # --- hierNet: hard 32-bit .C ceiling above ~1000 predictors ---------------
  if (want("hiernet") && has("hierNet")) {
    r = if (P > 1000)
          list(value=NULL, secs=NA_real_, peakMB=NA_real_, ok=FALSE,
               msg=sprintf("infeasible: p=%d exceeds the 32-bit .C limit", P))
        else timed({ setTimeLimit(elapsed = HN.BUDGET, transient = TRUE)
                     on.exit(setTimeLimit(elapsed = Inf, transient = TRUE))
                     A = do.call(cbind, Xtr); B = do.call(cbind, lapply(Xall, function(G) G[(nTrain+1):N,,drop=FALSE]))
                     fit = hierNet::hierNet.path(A, as.numeric(ytr), nlam = 5, trace = 0)
                     as.numeric(stats::predict(fit, newx = B)[,3]) })
    cc = if (r$ok) tryCatch(cor(yts, r$value), error=function(e) NA_real_) else NA_real_
    add.row("hierNet", r, cc); flush.rows()
  }
  flush.rows()
  rm(d, Xtr, Xall); gc(full = TRUE)
  cat(sprintf("  rep %d done\n", rep)); flush.console()
}

res = do.call(rbind, rows)
dir.create(file.path(RT, "result", "benchmark"), recursive = TRUE, showWarnings = FALSE)
write.csv(res, file.path(RT, "result", "benchmark",
          sprintf("bench_%s_N%d_P%d%s.csv", basis, N, nSNP, SUF)), row.names = FALSE)
print(res)
cat("[task ", id, "] DONE\n", sep = "")
