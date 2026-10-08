# ---------------------------------------------------------------------------
# Accuracy against sample size (Reviewer 2, comment 2).  Same style as
# sim_ngene_interaction_P20.pdf: one row of scenario panels, one line per
# method with +/- 1 SE ribbons, correlation only.
#
#   sim_nvary_nonlinear.pdf    1 x 5, the five link functions
#   sim_nvary_interaction.pdf  1 x 4, the four interaction mechanisms
#
# Fixed: 20 SNPs per gene, 10 regions, causal rate 0.2, h2 = 0.5.  x = sample
# size, on the same grid and with the same split rule as the runtime study
# (nTest = max(100, 0.2N)), so N = 1000 is the 800/200 split used throughout.
#
# N CEILING.  All nine architectures were run to N = 5000; the linear
# architecture was additionally run to N = 10,000 in an earlier pass.  The
# figure stops at NVARY_NMAX (default 5000) so every panel covers the same
# range; the linear result at 10,000 is reported in the text instead.
#
# Replicate counts are 100 at every level shown and are recorded per row in
# result/nvary_summary.csv rather than drawn, since they no longer vary.
# ---------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(ggpubr)})
theme_set(theme_grey(base_size = 19))
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
OUT = paste0(RT, "plot/")
NMAX = as.numeric(Sys.getenv("NVARY_NMAX","5000"))
pal = read.csv(paste0(RT,"scripts/method_palette.csv"), stringsAsFactors=FALSE)
PAL = setNames(pal$colour, pal$label); ORDER = pal$label

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

fs = list.files(paste0(RT,"result/sim_nvary"), pattern="\\.csv$", full.names=TRUE)
if (!length(fs)) stop("no sim_nvary results yet")
raw = do.call(rbind, lapply(fs, read.csv, stringsAsFactors=FALSE))
raw$Method = unname(canon[raw$method]); raw = raw[!is.na(raw$Method), ]
raw = raw[raw$N <= NMAX, ]
agg = do.call(rbind, lapply(split(raw, list(raw$scen, raw$N, raw$Method), drop=TRUE), function(z) {
  v = z$cor[is.finite(z$cor)]
  data.frame(scen=z$scen[1], N=z$N[1], Method=z$Method[1],
             mean=if(length(v)) mean(v) else NA_real_,
             se=if(length(v)>1) sd(v)/sqrt(length(v)) else NA_real_,
             usable=length(v), n=nrow(z)) }))
present = ORDER[ORDER %in% unique(agg$Method)]
agg$Method = factor(agg$Method, levels=present)
NS = sort(unique(agg$N))
YL = c(0, max(agg$mean + ifelse(is.na(agg$se),0,agg$se), na.rm=TRUE) * 1.04)
have = function(s) s %in% agg$scen

one = function(d) {
  ggplot(d, aes(N, mean, colour=Method, group=Method)) +
    geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Method), alpha=.15,
                colour=NA, show.legend=FALSE) +
    geom_line(size=.8) + geom_point(size=1.8) +
    scale_colour_manual(values=PAL[present], limits=present, breaks=present,
                        drop=FALSE, name="Method") +
    scale_fill_manual(values=PAL[present], limits=present, drop=FALSE, guide="none") +
    scale_x_log10(breaks=NS, labels=NS) + coord_cartesian(ylim=YL) +
    ggtitle(NICE[[d$scen[1]]]) +
    theme(plot.title=element_text(hjust=.5, size=16),
          axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x="Sample size", y="COR")
}
figure = function(scen, fname) {
  scen = scen[sapply(scen, have)]
  if (!length(scen)) { cat("skip", fname, "-- no results yet\n"); return(invisible()) }
  g = ggarrange(plotlist=lapply(scen, function(sc) one(agg[agg$scen==sc,])),
                nrow=1, ncol=length(scen), common.legend=TRUE, legend="bottom", align="hv")
  ggsave(paste0(OUT,fname), g, width=3.6*length(scen)+0.4, height=4.9,
         device="pdf", limitsize=FALSE)
  cat(sprintf("wrote %s (%d panels, N up to %d, ylim [%.2f, %.2f])\n",
              fname, length(scen), max(NS), YL[1], YL[2]))
  for (m in present) {
    z = agg[agg$Method==m & agg$scen %in% scen, ]
    lo = mean(z$mean[z$N==min(NS)], na.rm=TRUE); hi = mean(z$mean[z$N==max(NS)], na.rm=TRUE)
    cat(sprintf("   %-28s N=%d %.3f -> N=%d %.3f  (%+.3f)\n", m, min(NS), lo, max(NS), hi, hi-lo))
  }
}
figure(NONLIN, "sim_nvary_nonlinear.pdf")
figure(INTER,  "sim_nvary_interaction.pdf")
write.csv(agg, paste0(RT,"result/nvary_summary.csv"), row.names=FALSE)
cat("DONE\n")
