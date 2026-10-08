# ---------------------------------------------------------------------------
# Four focused simulation figures, all from the h2 = 0.5 arm, all 2 rows
# (MSE on top, COR below) x one column per scenario.
#
#   sim_h2_0.5_nonlinear_baselines.pdf    5 link functions  x {CF-ridge, GMM,
#                                         BLUP, LDpred2, hierNet}
#   sim_h2_0.5_interaction_baselines.pdf  4 interaction scenarios, same methods
#   sim_h2_0.5_nonlinear_minque.pdf       5 link functions  x the three MINQUE
#                                         variants
#   sim_h2_0.5_interaction_minque.pdf     4 interaction scenarios, same methods
#
# SNP DENSITY.  The requested 2 x ncol layout has rows = metric and columns =
# scenario, which leaves SNP density nowhere to go, so each figure is a single
# density.  It is 20 SNPs/gene because hierNet exists only at that density
# (above it, p exceeds R's 32-bit .C index limit), and pooling densities would
# put a 20-only hierNet box beside all-density boxes for everything else.
# Override with PANEL_PSNP=100 (etc.) for the two MINQUE figures, which have no
# such restriction.
#
# Colours come from plot/method_palette.csv, the same file the two large
# figures use, so a method is the same colour in every figure in this folder.
# y limits are per metric ROW -- shared across that row's panels so the columns
# are comparable -- and applied with coord_cartesian(), which zooms rather than
# dropping points the way ylim() would.
# ---------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(ggpubr)})
# Font sizes are set against the reduction these figures take on the page, not
# against how the PDF looks on its own: a ~18in figure placed at a 6.5in
# \textwidth shrinks by ~2.8x, so 19pt here lands near 7pt in the article.
# Column width is per scenario set, because the interaction labels are long.
theme_set(theme_grey(base_size = 19))
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
OUT = paste0(RT, "plot/")
ARM = "_h2_0.5"
PSNP = as.integer(Sys.getenv("PANEL_PSNP", "20"))

pal = read.csv(paste0(RT, "scripts/method_palette.csv"), stringsAsFactors = FALSE)
PAL = setNames(pal$colour, pal$label); ORDER = pal$label

canon = c("CLOSEFORM-MINQUE1-RIDGE"="CLOSEFORM-MINQUE-RIDGE",
  "NUMERICAL-MINQUE-LASSO"="NUMERICAL-MINQUE-LASSO",
  "NUMERICAL-MINQUE-ELASTICNET"="NUMERICAL-MINQUE-ELASTICNET",
  "GMM"="GMM", "BLUP"="BLUP", "LDpred2"="LDpred2", "hierNet"="hierNet")

CAUSAL = c("0.2","0.4","0.6","0.8","1")
NONLIN = c("linear","cosh","square","hyperbola","rickercurve")
INTER  = c("within.region.threshold.TRUE","within.region.threshold.FALSE",
           "outside.region.threshold.TRUE","outside.region.threshold.FALSE")
NICE = c(linear="Linear", cosh="Cosh", square="Square", hyperbola="Hyperbola",
  rickercurve="Rickercurve", within.region.threshold.TRUE="Threshold within Genes",
  within.region.threshold.FALSE="Multiplicative within Genes",
  outside.region.threshold.TRUE="Threshold between Genes",
  outside.region.threshold.FALSE="Multiplicative between Genes")

G.BASE   = c("CLOSEFORM-MINQUE-RIDGE","GMM","BLUP","LDpred2","hierNet")
G.MINQUE = c("CLOSEFORM-MINQUE-RIDGE","NUMERICAL-MINQUE-LASSO","NUMERICAL-MINQUE-ELASTICNET")

VALID = list()

cell = function(sf, cc) {
  f = sprintf("%sresult/sim%s/eval_%s_%s_%s.csv", RT, ARM, sf, PSNP, cc)
  if (!file.exists(f)) return(NULL)
  d = read.csv(f, stringsAsFactors = FALSE)[, c("rep","method","mse","cor")]
  for (t in c("ldpred2","hiernet100")) {
    g = sprintf("%sresult/sim%s_extra/eval_%s_%s_%s_%s.csv", RT, ARM, sf, PSNP, cc, t)
    if (file.exists(g))
      d = rbind(d, read.csv(g, stringsAsFactors = FALSE)[, c("rep","method","mse","cor")]) }
  d$causal = cc; d
}

# whisker span of every method in a set of panels, so nothing is drawn outside
# the limits (outliers are suppressed, so whiskers are what must be visible)
wrange = function(v, m) {
  r = do.call(rbind, lapply(split(v, m), function(x) {
        x = x[is.finite(x)]; if (length(x) < 2) return(c(NA, NA))
        st = grDevices::boxplot.stats(x)$stats; c(st[1], st[5]) }))
  lo = min(r[,1], na.rm = TRUE); hi = max(r[,2], na.rm = TRUE)
  lo = min(0, lo); c(lo, hi + 0.04 * max(hi - lo, .05))
}

panel.figure = function(scen, methods, fname, ttl, colw = 3.6) {
  dat = lapply(scen, function(s) {
    d = do.call(rbind, lapply(CAUSAL, function(cc) cell(s, cc)))
    if (is.null(d)) return(NULL)
    d$Method = unname(canon[d$method])
    d = d[!is.na(d$Method) & d$Method %in% methods, ]
    d$Causal = factor(d$causal, levels = CAUSAL); d })
  names(dat) = scen
  keep = !sapply(dat, is.null)
  if (!any(keep)) { cat("skip", fname, "-- no data\n"); return(invisible()) }
  dat = dat[keep]; scen = scen[keep]

  all = do.call(rbind, dat)
  present = ORDER[ORDER %in% unique(all$Method)]
  miss = setdiff(methods, present)
  if (length(miss)) cat("  NOTE", fname, "-- absent at", PSNP, "SNPs/gene:",
                        paste(miss, collapse = ", "), "\n")
  ylim.mse = wrange(all$mse, all$Method)
  ylim.cor = wrange(all$cor, all$Method)
  cat(sprintf("  %s: %d methods, %d scenarios, MSE [%.2f, %.2f], COR [%.2f, %.2f]\n",
      fname, length(present), length(scen), ylim.mse[1], ylim.mse[2],
      ylim.cor[1], ylim.cor[2]))

  one = function(d, metric, yl, title, xlab) {
    d$Method = factor(d$Method, levels = present)
    d$value  = d[[metric]]
    for (m in present) VALID[[length(VALID)+1]] <<- data.frame(figure = fname,
      panel = title, metric = toupper(metric), method = m,
      n_total = sum(d$Method == m),
      n_usable = sum(is.finite(d$value[d$Method == m])))
    ggplot(d, aes(Causal, value, fill = Method)) +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = yl) +
      scale_fill_manual(values = PAL[present], limits = present, breaks = present,
                        drop = FALSE, name = "Method") +
      (if (nzchar(title)) ggtitle(title) else ggtitle(NULL)) +
      theme(plot.title = element_text(hjust = .5, size = 17)) +
      labs(x = xlab, y = toupper(metric)) }

  panels = c(lapply(scen, function(s) one(dat[[s]], "mse", ylim.mse, NICE[[s]], NULL)),
             lapply(scen, function(s) one(dat[[s]], "cor", ylim.cor, "",  "Causal Rate")))
  # Legend along the bottom: at this font size the method names are too wide
  # for a right-hand legend without eating the panels.
  g = ggarrange(plotlist = panels, nrow = 2, ncol = length(scen),
                common.legend = TRUE, legend = "bottom", align = "hv")
  # No figure title: every descriptive detail lives in the LaTeX caption.
  # The per-column scenario labels stay -- they identify the columns.
  ggsave(paste0(OUT, fname), g, width = colw * length(scen) + 0.4, height = 6.8,
         device = "pdf", limitsize = FALSE)
  cat("wrote", fname, "\n"); flush.console()
}

panel.figure(NONLIN, G.BASE,   sprintf("sim_h2_0.5_nonlinear_baselines_P%d.pdf", PSNP),
             "Nonlinear architectures: PKNN vs external baselines", colw = 3.6)
panel.figure(INTER,  G.BASE,   sprintf("sim_h2_0.5_interaction_baselines_P%d.pdf", PSNP),
             "Interaction architectures: PKNN vs external baselines", colw = 4.3)
panel.figure(NONLIN, G.MINQUE, sprintf("sim_h2_0.5_nonlinear_minque_P%d.pdf", PSNP),
             "Nonlinear architectures: MINQUE penalty variants", colw = 3.6)
panel.figure(INTER,  G.MINQUE, sprintf("sim_h2_0.5_interaction_minque_P%d.pdf", PSNP),
             "Interaction architectures: MINQUE penalty variants", colw = 4.3)

write.csv(do.call(rbind, VALID), sprintf("%sresult/panel_validity_P%d.csv", RT, PSNP),
          row.names = FALSE)
cat("DONE\n")
