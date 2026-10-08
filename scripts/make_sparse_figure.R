# ---------------------------------------------------------------------------
# Within-region sparsity figures (Reviewer 2, comment 1).  Same style as
# sim_ngene_interaction_P20.pdf: one row of scenario panels, x = the swept
# quantity, one line per method with +/- 1 SE ribbons, correlation only.
#
#   sim_sparse_nonlinear.pdf    1 x 5, the five link functions
#   sim_sparse_interaction.pdf  1 x 4, the four interaction mechanisms
#
# Fixed: h2 = 0.5, 20 SNPs per gene, 10 regions.  x = the PROPORTION of SNPs
# inside each causal region that carry an effect: 5, 10, 15 and 20 of the 20
# SNPs, i.e. 0.25, 0.50, 0.75 and 1.00.  At 1.00 every SNP in a causal region is
# causal, which reproduces the design used elsewhere in the paper and is
# therefore the reference level.
#
# CAUSAL RATE.  The data cover rates 0.2-1.0, but accuracy varies strongly with
# that rate, so pooling would give ribbons dominated by it rather than by
# replicate noise.  The figure therefore shows one rate, 0.2 by default (two of
# ten regions causal, matching the sample-size study), set by SPARSE_CAUSAL.
# ---------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(ggpubr)})
theme_set(theme_grey(base_size = 19))
RT  = "/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/"
OUT = paste0(RT, "plot/")
pal = read.csv(paste0(RT,"scripts/method_palette.csv"), stringsAsFactors=FALSE)
PAL = setNames(pal$colour, pal$label); ORDER = pal$label

NCAUSAL = c(5,10,15,20)
CAUSAL  = Sys.getenv("SPARSE_CAUSAL","0.2")
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

cell = function(sc, nc) {
  d = list()
  f = sprintf("%sresult/sim_sparse/eval_%s_%d_%s.csv", RT, sc, nc, CAUSAL)
  if (file.exists(f)) d[[1]] = read.csv(f, stringsAsFactors=FALSE)[,c("rep","method","cor")]
  for (t in c("ldpred2","hiernet")) {
    g = sprintf("%sresult/sim_sparse_extra/eval_%s_%d_%s_%s.csv", RT, sc, nc, CAUSAL, t)
    if (file.exists(g)) d[[length(d)+1]] = read.csv(g, stringsAsFactors=FALSE)[,c("rep","method","cor")]
  }
  if (!length(d)) return(NULL)
  z = do.call(rbind, d); z$scen = sc; z$x = nc / 20; z   # proportion of the 20 SNPs
}
have = function(s) any(file.exists(sprintf("%sresult/sim_sparse/eval_%s_%d_%s.csv", RT, s, NCAUSAL, CAUSAL)))
SC = c(NONLIN, INTER)[sapply(c(NONLIN, INTER), have)]
if (!length(SC)) stop("no sim_sparse results for causal rate ", CAUSAL)

raw = do.call(rbind, lapply(SC, function(sc) do.call(rbind, lapply(NCAUSAL, function(nc) cell(sc,nc)))))
raw$Method = unname(canon[raw$method]); raw = raw[!is.na(raw$Method), ]
agg = do.call(rbind, lapply(split(raw, list(raw$scen, raw$x, raw$Method), drop=TRUE), function(z) {
  v = z$cor[is.finite(z$cor)]
  data.frame(scen=z$scen[1], x=z$x[1], Method=z$Method[1],
             mean=if(length(v)) mean(v) else NA_real_,
             se=if(length(v)>1) sd(v)/sqrt(length(v)) else NA_real_,
             usable=length(v)/nrow(z), n=nrow(z)) }))
present = ORDER[ORDER %in% unique(agg$Method)]
agg$Method = factor(agg$Method, levels=present)
YL = c(0, max(agg$mean + ifelse(is.na(agg$se),0,agg$se), na.rm=TRUE) * 1.04)

one = function(d) {
  ggplot(d, aes(x, mean, colour=Method, group=Method)) +
    geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Method), alpha=.15,
                colour=NA, show.legend=FALSE) +
    geom_line(size=.8) + geom_point(size=1.8) +
    scale_colour_manual(values=PAL[present], limits=present, breaks=present,
                        drop=FALSE, name="Method") +
    scale_fill_manual(values=PAL[present], limits=present, drop=FALSE, guide="none") +
    scale_x_continuous(breaks=NCAUSAL/20, labels=sprintf("%.2f", NCAUSAL/20)) +
    coord_cartesian(ylim=YL) +
    ggtitle(NICE[[d$scen[1]]]) +
    theme(plot.title=element_text(hjust=.5, size=16)) +
    # "Causal SNP proportion" (194pt of Helvetica at 19pt) fits a 3.6in
    # column; the full wording "within each causal region" was 306pt and
    # overran the panel.  The qualifier belongs in the caption.
    labs(x="Causal SNP proportion", y="COR")
}
figure = function(scen, fname) {
  scen = scen[scen %in% SC]
  if (!length(scen)) { cat("skip", fname, "-- no results\n"); return(invisible()) }
  g = ggarrange(plotlist=lapply(scen, function(sc) one(agg[agg$scen==sc,])),
                nrow=1, ncol=length(scen), common.legend=TRUE, legend="bottom", align="hv")
  ggsave(paste0(OUT,fname), g, width=3.6*length(scen)+0.4, height=4.9,
         device="pdf", limitsize=FALSE)
  cat(sprintf("wrote %s (%d panels, causal rate %s, ylim [%.2f, %.2f])\n",
              fname, length(scen), CAUSAL, YL[1], YL[2]))
  for (m in present) {
    z = agg[agg$Method==m & agg$scen %in% scen, ]
    cat(sprintf("   %-28s prop 0.25 %.3f -> 1.00 %.3f  (usable %.0f%%)\n", m,
        mean(z$mean[abs(z$x-0.25)<1e-9], na.rm=TRUE),
        mean(z$mean[abs(z$x-1.00)<1e-9], na.rm=TRUE), 100*mean(z$usable)))
  }
}
figure(NONLIN, "sim_sparse_nonlinear.pdf")
figure(INTER,  "sim_sparse_interaction.pdf")
write.csv(agg, paste0(RT,"result/sparse_summary.csv"), row.names=FALSE)
cat("DONE\n")
