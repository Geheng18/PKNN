# ---------------------------------------------------------------------------
# ADNI result tables (Experiments 1 and 2), computed from the stored per-split
# results rather than transcribed.  Emits CSV, a plain-text/markdown copy for
# pasting, and a rendered PDF table.
#   value  = mean correlation (SD) over the same 100 splits
#   dagger = splits where the genetic fit collapsed, so the score is the
#            covariate prediction alone
#   approx = statistically tied with the best method in that column (paired t,
#            p > 0.05, same splits)
# CF-MINQUE1-ridge uses lambda chosen by 5-fold CV inside each training set,
# on unnormalised kernels.
# ---------------------------------------------------------------------------
suppressMessages({library(ggpubr)})
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
OUT = paste0(RT, "plot/"); dir.create(OUT, showWarnings = FALSE)
PH = c("AV45","Entorhinal","FDG","Hippocampus")

gather = function(expdir, cvdir) {
  do.call(rbind, lapply(PH, function(ph) {
    m  = read.csv(sprintf("%srealdata/%s/result/%s.csv", RT, expdir, ph), stringsAsFactors = FALSE)
    m  = m[m$method != "LDpred2", c("phenotype","rep","method","cor")]
    lp = read.csv(sprintf("%srealdata/%s/result/%s_ldpred2.csv", RT, expdir, ph), stringsAsFactors = FALSE)
    cv = read.csv(sprintf("%srealdata/lambda_cv/%s/result/%s.csv", RT, cvdir, ph), stringsAsFactors = FALSE)
    cv = cv[cv$kernels == "raw" & cv$lambda.rule == "CV", ]
    rbind(m, lp[, c("phenotype","rep","method","cor")],
          data.frame(phenotype = ph, rep = cv$rep, method = "CF-MINQUE1-RIDGE-CV", cor = cv$cor_final))
  }))
}

build = function(expdir, cvdir, rows, labels, file.stem, caption) {
  d = gather(expdir, cvdir)
  tab = matrix("", length(rows), length(PH), dimnames = list(labels, PH))
  collapsed = matrix(0L, length(rows), length(PH), dimnames = list(labels, PH))
  for (j in seq_along(PH)) {
    ph = PH[j]; sub = d[d$phenotype == ph, ]
    cov = sub[sub$method == "COVARIATE-ONLY", c("rep","cor")]
    mu = sapply(rows, function(m) mean(sub$cor[sub$method == m], na.rm = TRUE))
    best = rows[which.max(mu)]
    bv = sub[sub$method == best, c("rep","cor")]
    for (i in seq_along(rows)) {
      x = sub[sub$method == rows[i], c("rep","cor")]
      me = mean(x$cor, na.rm = TRUE); sd. = sd(x$cor, na.rm = TRUE)
      cl = if (rows[i] == "COVARIATE-ONLY") 0L else {
        z = merge(x, cov, by = "rep"); sum(abs(z$cor.x - z$cor.y) < 1e-9, na.rm = TRUE) }
      collapsed[i, j] = cl
      tie = FALSE
      if (rows[i] != best) {
        z = merge(x, bv, by = "rep"); df = z$cor.x - z$cor.y; ok = is.finite(df)
        tie = sum(ok) > 2 && sd(df[ok]) > 0 && t.test(df[ok])$p.value > 0.05
      }
      tab[i, j] = sprintf("%.3f (%.3f)%s%s", me, sd.,
                          if (cl > 0) "†" else "", if (tie) " ≈" else "")
    }
  }
  out = data.frame(method = labels, tab, check.names = FALSE)
  write.csv(out, paste0(OUT, file.stem, ".csv"), row.names = FALSE)
  cc = data.frame(method = labels, collapsed, check.names = FALSE)
  write.csv(cc, paste0(OUT, file.stem, "_collapsed.csv"), row.names = FALSE)
  # plain text for pasting
  con = file(paste0(OUT, file.stem, ".txt"), "w")
  writeLines(caption, con)
  writeLines(paste(c("method", PH), collapse = "\t"), con)
  for (i in seq_along(labels)) writeLines(paste(c(labels[i], tab[i, ]), collapse = "\t"), con)
  writeLines("", con)
  writeLines(paste0("Collapsed splits out of 100 (", paste(PH, collapse = " / "), "):"), con)
  for (i in seq_along(labels)) if (any(collapsed[i, ] > 0))
    writeLines(paste0("  ", labels[i], " ", paste(collapsed[i, ], collapse = " / ")), con)
  close(con)
  g = ggtexttable(out, rows = NULL, theme = ttheme("light", base_size = 9))
  ggsave(paste0(OUT, file.stem, ".pdf"), g, width = 11, height = 3.2, device = "pdf")
  cat("wrote", file.stem, ".csv/.txt/.pdf\n"); flush.console()
  out
}

R1 = c("COVARIATE-ONLY","CF-MINQUE1-RIDGE-CV","NUM-LASSO","NUM-ENET","GMM","BLUP","LDpred2")
L1 = c("Covariates only","CF-MINQUE1-ridge","NUM-lasso","NUM-enet",
       "GMM (CV lambda, as published)","BLUP (candidate GRM)","LDpred2 (10,010 candidate SNPs)")
L2 = c("Covariates only","CF-MINQUE1-ridge (full interaction)","NUM-lasso","NUM-enet",
       "GMM (CV lambda, as published)","BLUP (whole-genome GRM)","LDpred2 (649,072 pruned SNPs)")
t1 = build("exp1_candidate",  "exp1_candidate",  R1, L1, "table_exp1_candidate",
           "Experiment 1: 11 candidate genes")
t2 = build("exp2_genomewide", "exp2_genomewide", R1, L2, "table_exp2_genomewide",
           "Experiment 2: candidate genes plus the rest of the genome")
print(t1); print(t2)
cat("DONE\n")
