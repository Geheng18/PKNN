# ---------------------------------------------------------------------------
# "How many candidate genes?" figures (AE comment 1.4).  Two figures:
#
#   sim_ngene_nonlinear_P20.pdf    1 x 5, the five link functions
#   sim_ngene_interaction_P20.pdf  1 x 4, the four interaction mechanisms
#
# h2 = 0.5, 20 SNPs per gene, 2 CAUSAL regions throughout; the number of input
# regions grows from 2 (the oracle, only the causal genes supplied) to 20, the
# ceiling imposed by one region per chromosome on the UK Biobank panel.
#
# Lines rather than boxplots: the message is the trend in the number of genes
# supplied, and 8 region counts x 7 methods would be 56 boxes per panel.
# Ribbons are +/- 1 standard error of the mean over replicates.  Correlation
# only; retention against the R = 2 oracle is a rescaling of the same numbers
# and is kept in result/ngene_summary.csv for the text.
#
# The y range is computed over BOTH arms and shared, so the two figures can be
# read against each other.
#
# All lines solid, with a single Method legend.  The usable-fit rate is not
# encoded in the figure: it is reported per method in result/ngene_summary.csv
# and belongs in the caption.  GMM is usable in roughly a third of replicates,
# so its apparent improvement with more regions is not robustness.
#
# Colours are the shared palette, so a method keeps its colour across every
# figure in plot/.
# ---------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(ggpubr)})
# ~18in wide and reducing ~2.8x at a 6.5in \textwidth, so 19pt lands near 7pt.
theme_set(theme_grey(base_size = 19))
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
OUT = paste0(RT, "plot/")
pal = read.csv(paste0(RT,"scripts/method_palette.csv"), stringsAsFactors=FALSE)
PAL = setNames(pal$colour, pal$label); ORDER = pal$label

REGION = c(2,4,6,8,10,12,16,20)
NONLIN = c("linear","cosh","square","hyperbola","rickercurve")
INTER  = c("within.region.threshold.TRUE","within.region.threshold.FALSE",
           "outside.region.threshold.TRUE","outside.region.threshold.FALSE")
NICE   = c(linear="Linear", cosh="Cosh", square="Square",
           hyperbola="Hyperbola", rickercurve="Rickercurve",
           within.region.threshold.TRUE="Threshold within Genes",
           within.region.threshold.FALSE="Multiplicative within Genes",
           outside.region.threshold.TRUE="Threshold between Genes",
           outside.region.threshold.FALSE="Multiplicative between Genes")
canon = c("CLOSEFORM-MINQUE1-RIDGE"="CLOSEFORM-MINQUE-RIDGE",
          "NUMERICAL-MINQUE-LASSO"="NUMERICAL-MINQUE-LASSO",
          "NUMERICAL-MINQUE-ELASTICNET"="NUMERICAL-MINQUE-ELASTICNET",
          "GMM"="GMM","BLUP"="BLUP","LDpred2"="LDpred2","hierNet"="hierNet")

cell = function(sc, r) {
  d = list()
  f = sprintf("%sresult/sim_ngene/eval_%s_%d.csv", RT, sc, r)
  if (file.exists(f)) d[[1]] = read.csv(f, stringsAsFactors=FALSE)[,c("rep","method","cor")]
  for (t in c("ldpred2","hiernet")) {
    g = sprintf("%sresult/sim_ngene_extra/eval_%s_%d_%s.csv", RT, sc, r, t)
    if (file.exists(g)) d[[length(d)+1]] = read.csv(g, stringsAsFactors=FALSE)[,c("rep","method","cor")]
  }
  if (!length(d)) return(NULL)
  z = do.call(rbind, d); z$scen = sc; z$R = r; z
}
have = function(s) file.exists(sprintf("%sresult/sim_ngene/eval_%s_2.csv", RT, s))
SC   = c(NONLIN, INTER)[sapply(c(NONLIN, INTER), have)]
if (!length(SC)) stop("no sim_ngene results found")

raw = do.call(rbind, lapply(SC, function(sc) do.call(rbind, lapply(REGION, function(r) cell(sc,r)))))
raw$Method = unname(canon[raw$method]); raw = raw[!is.na(raw$Method), ]

agg = do.call(rbind, lapply(split(raw, list(raw$scen, raw$R, raw$Method), drop=TRUE), function(z) {
  v = z$cor[is.finite(z$cor)]
  data.frame(scen=z$scen[1], R=z$R[1], Method=z$Method[1],
             mean=if(length(v)) mean(v) else NA_real_,
             se  =if(length(v)>1) sd(v)/sqrt(length(v)) else NA_real_,
             usable=length(v)/nrow(z), n=nrow(z)) }))
# retention against each method's own R = 2 value, kept for the text
agg$ret = NA_real_
for (sc in SC) for (m in unique(agg$Method)) {
  i = agg$scen==sc & agg$Method==m
  b = agg$mean[i & agg$R==2]
  if (length(b)==1 && is.finite(b) && b > 0) agg$ret[i] = 100*agg$mean[i]/b
}
present = ORDER[ORDER %in% unique(agg$Method)]
agg$Method = factor(agg$Method, levels=present)
# one y range over both arms
YL = c(0, max(agg$mean + ifelse(is.na(agg$se),0,agg$se), na.rm=TRUE) * 1.04)

one = function(d) {
  ggplot(d, aes(R, mean, colour=Method, group=Method)) +
    geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Method), alpha=.15,
                colour=NA, show.legend=FALSE) +
    geom_line(size=.8) + geom_point(size=1.8) +
    scale_colour_manual(values=PAL[present], limits=present, breaks=present,
                        drop=FALSE, name="Method") +
    scale_fill_manual(values=PAL[present], limits=present, drop=FALSE, guide="none") +
    scale_x_continuous(breaks=REGION) + coord_cartesian(ylim=YL) +
    ggtitle(NICE[[d$scen[1]]]) +
    # 16pt, not 17: the widest title ("Multiplicative between Genes") is
    # 221pt of Helvetica at 17pt and ran 4pt past the right edge of the
    # four-panel interaction figure.
    theme(plot.title=element_text(hjust=.5, size=16)) +
    labs(x="Number of gene regions", y="COR")
}
figure = function(scen, fname) {
  scen = scen[scen %in% SC]
  if (!length(scen)) { cat("skip", fname, "-- no results yet\n"); return(invisible()) }
  panels = lapply(scen, function(sc) one(agg[agg$scen==sc,]))
  # The axis label is each panel's own x title, which ggplot always draws
  # between the panel and a bottom legend.  Two attempts to place it by
  # composing rows -- annotate_figure(bottom=) and an explicit 3-row
  # ggarrange -- both rendered it UNDERNEATH the legend.  Shortening the label
  # to "Number of gene regions" (204pt of Helvetica at 19pt) makes it fit a
  # single ~220pt column, which the longer wording did not.
  g = ggarrange(plotlist=panels, nrow=1, ncol=length(scen),
                common.legend=TRUE, legend="bottom", align="hv")
  ggsave(paste0(OUT,fname), g, width=3.6*length(scen)+0.4, height=4.9,
         device="pdf", limitsize=FALSE)
  cat(sprintf("wrote %s (%d panels, ylim [%.2f, %.2f])\n", fname, length(scen), YL[1], YL[2]))
  for (m in present) {
    z = agg[agg$Method==m & agg$scen %in% scen, ]
    cat(sprintf("   %-28s R=2 %.3f -> R=20 %.3f  (%.0f%% retained, usable %.0f%%)\n", m,
        mean(z$mean[z$R==2], na.rm=TRUE), mean(z$mean[z$R==20], na.rm=TRUE),
        mean(z$ret[z$R==20], na.rm=TRUE), 100*mean(z$usable)))
  }
}
figure(NONLIN, "sim_ngene_nonlinear_P20.pdf")
figure(INTER,  "sim_ngene_interaction_P20.pdf")
write.csv(agg, paste0(RT,"result/ngene_summary.csv"), row.names=FALSE)
cat("DONE\n")
