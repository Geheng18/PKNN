# Reviewer-requested comparators.
#   glinternet  - pairwise-interaction group lasso with strong hierarchy (CRAN)
#   bigsnpr     - LDpred2 (CRAN)
#   hierNet     - hierarchical interaction lasso, strong/weak heredity (CRAN).
#                 Maintained stand-in for SHIM, which was never on CRAN.
#   shim        - Choi/Li/Zhu strong-heredity model; GitHub only, last commit 2017.
#                 Attempted last, failure is tolerated.
LIB = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/Rlib"
dir.create(LIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)

for (p in c("remotes","glinternet","hierNet","bigsnpr")) {
  cat("\n========== installing", p, "==========\n"); flush.console()
  try(install.packages(p, lib = LIB))
  cat(p, ":", if (requireNamespace(p, quietly=TRUE)) "OK" else "FAILED", "\n"); flush.console()
}

cat("\n========== attempting shim from GitHub (unmaintained since 2017) ==========\n")
try(remotes::install_github("sahirbhatnagar/shim", lib = LIB, upgrade = "never"))

cat("\n\n=============== SUMMARY ===============\n")
for (p in c("glinternet","hierNet","bigsnpr","shim")) {
  ok = requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-12s %s %s\n", p, if (ok) "OK     " else "MISSING",
              if (ok) as.character(packageVersion(p)) else ""))
}
