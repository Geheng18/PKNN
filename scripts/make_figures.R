# ---------------------------------------------------------------------------
# Final figure set. Five figures, correlation only.
#   sim_h2_0.2_correlation.pdf   one figure, every scenario x SNP density
#   sim_h2_0.5_correlation.pdf   same for the main heritability arm
#   realdata_exp1.pdf            ADNI, 11 candidate genes
#   realdata_exp2.pdf            ADNI, candidate genes + rest of genome
#   runtime_memory.pdf           time and peak memory vs sample size, one figure
#
# One fixed colour per method (plot/method_palette.csv), shared across every
# figure.  Legends drop methods absent from the figure.
# Per-panel y limits are set from the 2nd-98th percentile of the values in that
# panel, padded, via coord_cartesian -- it zooms, so no data is removed from the
# boxplot statistics (ylim() would delete out-of-range points).
# Interaction scenarios are picked up automatically once their results exist.
# ---------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(ggpubr)})
# Font sizes are set against the reduction each figure takes on the page rather
# than against the PDF on its own; each ggsave below notes its own budget.
theme_set(theme_grey(base_size = 20))
RT = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
OUT = paste0(RT, "plot/"); dir.create(OUT, showWarnings = FALSE)
pal = read.csv(paste0(RT, "scripts/method_palette.csv"), stringsAsFactors = FALSE)
PAL = setNames(pal$colour, pal$label); ORDER = pal$label
canon = c("CLOSEFORM-MINQUE1-RIDGE"="CLOSEFORM-MINQUE-RIDGE","CF-MINQUE1-RIDGE-CV"="CLOSEFORM-MINQUE-RIDGE",
  "NUMERICAL-MINQUE-LASSO"="NUMERICAL-MINQUE-LASSO","NUM-LASSO"="NUMERICAL-MINQUE-LASSO",
  "NUMERICAL-MINQUE-ELASTICNET"="NUMERICAL-MINQUE-ELASTICNET","NUM-ENET"="NUMERICAL-MINQUE-ELASTICNET",
  "GMM"="GMM","BLUP"="BLUP","LDpred2"="LDpred2","glinternet"="glinternet","hierNet"="hierNet",
  "COVARIATE-ONLY"="Covariates only")
# glinternet dropped as a benchmark at the author's request.  Its results remain
# in result/sim*_extra/ and can be restored by adding it back to this vector.
SIM.METHODS  = c("CLOSEFORM-MINQUE-RIDGE","NUMERICAL-MINQUE-LASSO","NUMERICAL-MINQUE-ELASTICNET",
                 "GMM","BLUP","LDpred2","hierNet")
REAL.METHODS = c("Covariates only","CLOSEFORM-MINQUE-RIDGE","NUMERICAL-MINQUE-LASSO",
                 "NUMERICAL-MINQUE-ELASTICNET","GMM","BLUP","LDpred2")
CAUSAL = c("0.2","0.4","0.6","0.8","1"); PSNP = c(20,100,500,1000)
NONLIN = c("linear","cosh","square","hyperbola","rickercurve")
INTER  = c("within.region.threshold.TRUE","within.region.threshold.FALSE",
           "outside.region.threshold.TRUE","outside.region.threshold.FALSE")
NICE = c(linear="Linear", cosh="Cosh", square="Square", hyperbola="Hyperbola",
  rickercurve="Rickercurve", within.region.threshold.TRUE="Threshold within Genes",
  within.region.threshold.FALSE="Multiplicative within Genes",
  outside.region.threshold.TRUE="Threshold between Genes",
  outside.region.threshold.FALSE="Multiplicative between Genes")
VALID = list()
note = function(fig, panel, d) for (m in unique(d$Method))
  VALID[[length(VALID)+1]] <<- data.frame(figure=fig, panel=panel, method=m,
    n_total=sum(d$Method==m), n_usable=sum(is.finite(d$cor[d$Method==m])))
# One y range for the whole figure, taken from the boxplot whiskers of every
# method in every panel so nothing is drawn outside it, and anchored at 0 unless
# the data actually goes negative.  Outliers are not shown (outlier.shape = NA),
# so whisker extent is what has to be visible.
wrange = function(d) {
  r = do.call(rbind, lapply(split(d$cor, d$Method), function(x) {
        x = x[is.finite(x)]; if (length(x) < 2) return(c(NA, NA))
        st = grDevices::boxplot.stats(x)$stats; c(st[1], st[5]) }))
  c(min(r[,1], na.rm=TRUE), max(r[,2], na.rm=TRUE)) }
fig.ylim = function(lo, hi) { lo = min(0, lo); c(lo, hi + 0.04*max(hi-lo, .05)) }

sim.cell = function(arm, sf, p, cc) {
  f = sprintf("%sresult/sim%s/eval_%s_%s_%s.csv", RT, arm, sf, p, cc)
  if (!file.exists(f)) return(NULL)
  d = read.csv(f, stringsAsFactors=FALSE)[, c("rep","method","cor")]
  for (t in c("ldpred2","glinternet100","hiernet100")) {
    g = sprintf("%sresult/sim%s_extra/eval_%s_%s_%s_%s.csv", RT, arm, sf, p, cc, t)
    if (file.exists(g)) d = rbind(d, read.csv(g, stringsAsFactors=FALSE)[, c("rep","method","cor")]) }
  d$causal = cc; d
}
sim.figure = function(arm, fname) {
  scen = c(NONLIN, INTER)
  scen = scen[sapply(scen, function(s) file.exists(sprintf("%sresult/sim%s/eval_%s_20_0.2.csv", RT, arm, s)))]
  # pass 1: read every cell, so the legend and the y range can be set from the
  # whole figure rather than from whichever panel happens to be drawn first
  # Scenario-major: ggarrange fills row by row, so iterating scenario first
  # puts one architecture on each row and the four SNP densities across it.
  dat = list(); k = 0
  for (s in scen) for (p in PSNP) {
    d = do.call(rbind, lapply(CAUSAL, function(cc) sim.cell(arm, s, p, cc)))
    k = k + 1
    if (!is.null(d)) {
      d$Method = unname(canon[d$method]); d = d[!is.na(d$Method) & d$Method %in% SIM.METHODS, ]
      d$Causal = factor(d$causal, levels = CAUSAL)
    }
    dat[[k]] = list(d = d, s = s, p = p)
  }
  present = ORDER[ORDER %in% unique(unlist(lapply(dat, function(z) z$d$Method)))]
  allrows = do.call(rbind, lapply(dat, function(z) if (is.null(z$d)) NULL else z$d[, c("Method","cor")]))
  allrows$Method = factor(allrows$Method, levels = present)
  yl = do.call(fig.ylim, as.list(wrange(allrows)))
  cat(sprintf("  %s: %d methods in legend, ylim [%.2f, %.2f]\n", fname, length(present), yl[1], yl[2]))

  panels = lapply(dat, function(z) {
    if (is.null(z$d)) return(ggplot() + theme_void())
    d = z$d; d$Method = factor(d$Method, levels = present)
    note(fname, paste0(NICE[z$s], " P=", z$p), d)
    ggplot(d, aes(Causal, cor, fill=Method)) + geom_boxplot(outlier.shape=NA) +
      coord_cartesian(ylim = yl) +
      scale_fill_manual(values = PAL[present], limits = present, breaks = present,
                        drop = FALSE, name = "Method") +
      ggtitle(sprintf("%s  \u00b7  %d SNPs/gene", NICE[z$s], z$p)) +
      # 13pt in a 4.6in panel.  With only four columns the widest titles no
      # longer have a long neighbour to collide with, but keep headroom.
      theme(plot.title=element_text(hjust=.5, size=13)) + labs(x="Causal Rate", y="COR") })
  g = ggarrange(plotlist=panels, nrow=length(scen), ncol=length(PSNP),
                common.legend=TRUE, legend="right")
  # 9 rows x 4 columns is taller than it is wide, so the page constrains this
  # figure by height: ~22.6 x 31.5in placed at 0.9\textheight comes out about
  # 6.4in wide, which fits \textwidth on a portrait page.  The old 4 x 9
  # layout had to go sideways and still reduced harder.
  ggsave(paste0(OUT,fname), g, width=4.6*length(PSNP)+4.2, height=3.5*length(scen),
         device="pdf", limitsize=FALSE)
  cat("wrote", fname, "(", length(scen), "scenario rows x", length(PSNP), "density columns )\n"); flush.console()
}
sim.figure("",         "sim_h2_0.2_correlation.pdf")
sim.figure("_h2_0.5",  "sim_h2_0.5_correlation.pdf")

PH = c("AV45","Entorhinal","FDG","Hippocampus")
real.figure = function(expdir, cvdir, fname, title) {
  d = do.call(rbind, lapply(PH, function(ph) {
    m  = read.csv(sprintf("%srealdata/%s/result/%s.csv", RT, expdir, ph), stringsAsFactors=FALSE)
    m  = m[m$method != "LDpred2" & !grepl("^CF-", m$method), c("phenotype","rep","method","cor")]
    lp = read.csv(sprintf("%srealdata/%s/result/%s_ldpred2.csv", RT, expdir, ph), stringsAsFactors=FALSE)
    cv = read.csv(sprintf("%srealdata/lambda_cv/%s/result/%s.csv", RT, cvdir, ph), stringsAsFactors=FALSE)
    cv = cv[cv$kernels=="raw" & cv$lambda.rule=="CV", ]
    rbind(m, lp[,c("phenotype","rep","method","cor")],
          data.frame(phenotype=ph, rep=cv$rep, method="CF-MINQUE1-RIDGE-CV", cor=cv$cor_final)) }))
  d$Method = unname(canon[d$method]); d = d[!is.na(d$Method) & d$Method %in% REAL.METHODS, ]
  present = ORDER[ORDER %in% unique(d$Method)]
  d$Method = factor(d$Method, levels = present)
  d$Phenotype = factor(d$phenotype, levels=PH)
  for (ph in PH) note(fname, ph, d[d$phenotype==ph, ])
  yl = do.call(fig.ylim, as.list(wrange(d)))
  cat(sprintf("  %s: %d methods in legend, ylim [%.2f, %.2f]\n", fname, length(present), yl[1], yl[2]))
  g = ggplot(d, aes(Phenotype, cor, fill=Method)) + geom_boxplot(outlier.shape=NA) +
    coord_cartesian(ylim = yl) +
    scale_fill_manual(values=PAL[present], limits=present, breaks=present,
                      drop=FALSE, name="Method") +
    labs(y="COR") +  # no title: descriptive detail lives in the LaTeX caption
    theme(legend.position="bottom", legend.text=element_text(size=15),
          legend.title=element_text(size=16)) +
    guides(fill = guide_legend(nrow = 3))
  # 12 x 6.8in at a 6.5in \textwidth is a 1.85x reduction, so the 20pt base
  # lands near 11pt on the page.  Three legend rows, not two: with two, the
  # layout needs four columns and "Covariates only" ran 17pt off the page.
  ggsave(paste0(OUT,fname), g, width=12, height=6.8, device="pdf")
  cat("wrote", fname, "\n"); flush.console()
}
real.figure("exp1_candidate","exp1_candidate","realdata_exp1.pdf","ADNI: 11 candidate genes")
real.figure("exp2_genomewide","exp2_genomewide","realdata_exp2.pdf",
            "ADNI: candidate genes + rest of genome")

# Runtime and memory across every method, sample size and SNP density.
# Two metric rows (time, peak memory) x one column per method; the line COLOUR
# encodes SNP density.  Note this is the one figure where colour does not mean
# method -- method is the panel, because 7 methods x 6 densities is 42 lines and
# would be unreadable on a single pair of axes.  Log scales on both metrics: the
# values span roughly five orders of magnitude.
bm = do.call(rbind, lapply(list.files(paste0(RT,"result/benchmark"), pattern="\\.csv$",
                                      full.names=TRUE), read.csv))
bm = bm[bm$ok %in% c(TRUE,"TRUE"), ]
# Panel labels and panel order are the same names and the same order the
# simulation and real-data legends use (plot/method_palette.csv), so a method
# is called one thing across every figure in this folder.
bench.canon = c("CF-MINQUE1-RIDGE"="CLOSEFORM-MINQUE-RIDGE",
                "NUM-MINQUE-LASSO"="NUMERICAL-MINQUE-LASSO",
                "NUM-MINQUE-ENET" ="NUMERICAL-MINQUE-ELASTICNET",
                "GMM"="GMM", "BLUP"="BLUP", "LDpred2"="LDpred2", "hierNet"="hierNet")
bm$M = unname(bench.canon[bm$method]); bm = bm[!is.na(bm$M), ]
stopifnot(all(unique(bm$M) %in% ORDER))
bm$M = factor(bm$M, levels = ORDER[ORDER %in% unique(bm$M)])
agg = aggregate(cbind(secs, peakMB) ~ N + nSNP + M, bm, mean)
agg$SNPs = factor(agg$nSNP, levels = sort(unique(agg$nSNP)))
long = rbind(data.frame(agg[,c("N","M","SNPs")], metric="Time (s)",          value=agg$secs),
             data.frame(agg[,c("N","M","SNPs")], metric="Peak memory (GB)",  value=agg$peakMB/1024))
long$metric = factor(long$metric, levels=c("Time (s)","Peak memory (GB)"))
cat(sprintf("  runtime_memory.pdf: %d methods, %d SNP densities, %d sample sizes\n",
            nlevels(bm$M), nlevels(agg$SNPs), length(unique(agg$N))))
g = ggplot(long, aes(N, value, colour=SNPs, group=SNPs)) +
      geom_line() + geom_point(size=1) +
      # The canonical names are long; break them after "MINQUE" so a strip
      # label fits its panel instead of running into the neighbouring one.
      facet_grid(metric ~ M, scales="free_y",
                 labeller = labeller(M = function(x) sub("-MINQUE-", "-MINQUE\n", x))) +
      scale_x_log10() + scale_y_log10() +
      scale_colour_viridis_d(name="SNPs per gene", option="viridis", end=.92) +
      labs(x="Individual Size") +
      theme(strip.text=element_text(size=13), legend.position="right")
ggsave(paste0(OUT,"runtime_memory.pdf"), g, width=16, height=7.4, device="pdf")
cat("wrote runtime_memory.pdf\n")
write.csv(do.call(rbind, VALID), paste0(RT,"result/figure_validity.csv"), row.names=FALSE)
cat("DONE\n")
