# Code to reproduce figures, table, & key summary statistics
# Assumes the MC simulation code has been run, so that the
# .Rdata files with simulation results exist.
#
# Author: Steve Miller (steve.j.miller@colorado.edu)

# Load libraries
library(here)
library(ggplot2)
library(cowplot)
library(kableExtra)

# Helper function to summarize MC output
SummarizeMCRuns = function(mc.exp.list) {
  rbindlist(lapply(mc.exp.list,
                   FUN=function(mc.res) {
                     n = mc.res$n[1]
                     p = mc.res$p[1]
                     if(!is.null(mc.res$k)) {
                       k = mc.res$k[1]  
                     } else {
                       k=5
                     }
                     if(is.null(mc.res$het)) {
                       het = F
                     } else {
                       het = mc.res$het[1]
                     }
                     
                     alpha = mc.res$alpha[1]
                     pdl.bch.est = mean(mc.res$bch.pdl.est)
                     pdl.est = mean(mc.res$pdl.est)
                     pdl.bias.perc = (pdl.est-alpha)/alpha
                     offset.est = mean(mc.res$pdl.offset.est)
                     bias.perc = (offset.est-alpha)/alpha
                     bias.perc.of.sd = (offset.est-alpha)/sd(mc.res$pdl.offset.est)
                     coverage = mean(mc.res$pdl.offset.cover)
                     coverage.auxadj = mean(mc.res$pdl.offset.cover.auxadj)
                     coverage.auxadjhc = mean(mc.res$pdl.offset.cover.auxadjhc)
                     coverage.cjn = mean(mc.res$pdl.offset.cover.cjn)
                     coverage.bch.pdl = mean(mc.res$bch.pdl.cover)
                     sd.bch.pdl.est = sd(mc.res$bch.pdl.est)
                     sd.offset.est = sd(mc.res$pdl.offset.est)
                     sd.pdl.est = sd(mc.res$pdl.est)
                     efficiency = sd(mc.res$pdl.offset.est)/sd(mc.res$allcontrol.est)
                     se = mean(mc.res$pdl.offset.se)
                     se.auxadj = mean(mc.res$pdl.offset.se.auxadj)
                     se.auxadjhc = mean(mc.res$pdl.offset.se.auxadjhc)
                     se.cjn = mean(mc.res$pdl.offset.se.cjn)
                     se.allcontrol = mean(mc.res$allcontrol.se)
                     se.bch.pdl = mean(mc.res$bch.pdl.se)
                     coverage.allcontrol = mean(mc.res$allcontrol.cover)
                     coverage.allcontrol.cjn = mean(mc.res$allcontrol.cover.cjn)
                     data.table(alpha, n, p, k, het, pdl.bch.est, pdl.est, pdl.bias.perc, offset.est, bias.perc, 
                                coverage, coverage.auxadj, coverage.auxadjhc, coverage.cjn, coverage.allcontrol, coverage.allcontrol.cjn, coverage.bch.pdl, 
                                se,  se.auxadj, se.auxadjhc, se.cjn, se.allcontrol, se.bch.pdl, 
                                sd.bch.pdl.est, sd.pdl.est, sd.offset.est, efficiency)
                   }))
}


# Load all simulation results
load(here("results","simpleres.Rdata"))
load(here("results","rsqres.Rdata"))
load(here("results","alphares.Rdata"))
load(here("results","pres.Rdata"))
load(here("results","kres.Rdata"))

#####################
# Figure 1 and associated summary statistics
#####################
# Make Figure 1
mc.plots = lapply(list(mc.res._0.5, mc.res.0, mc.res.0.5), FUN=function(mc.res) {
  alpha = mc.res$alpha[1]
  ggplot(data=mc.res, 
         mapping=aes(x=pdl.est)) + 
    geom_histogram(fill="#cccccc") +
    geom_vline(aes(color="True parameter", xintercept = alpha)) +
    geom_vline(aes(color="Mean estimate", xintercept = mean(mc.res$pdl.est))) +
    xlab(expression(hat(alpha))) + 
    ylab("# Monte Carlo samples") + 
    scale_color_manual(name="", 
                       breaks=c("True parameter", "Mean estimate"), 
                       values = c("blue", "red")) +
    theme_minimal() + 
    theme(legend.position = "bottom") + 
    ggtitle(bquote(alpha == .(alpha)))
})
leg = get_legend(mc.plots[[1]])
plot_grid(plot_grid(plotlist = lapply(mc.plots, FUN=function(p) { p + theme(legend.position="none")}), 
                    ncol = 3),
          plot_grid(leg),
          nrow=2,
          rel_heights = c(0.9, 0.1)) +
  theme(plot.margin = unit(c(0.5,0.1,0.1,0.1), "cm"))
ggsave(here("figures","simpleex.png"), width=6.5, height=2.5)

# Point estimates for modified estimator:
mean(mc.res._0.5$pdl.offset.est)
mean(mc.res.0$pdl.offset.est)
mean(mc.res.0.5$pdl.offset.est)

# Document that in all cases, all 5 relevant controls are selected, so this
# is not about underselection
summary(mc.res._0.5[,n.corr.pdl.controls])
summary(mc.res.0[,n.corr.pdl.controls])
summary(mc.res.0.5[,n.corr.pdl.controls])

# how much overselection is happening?
median(mc.res._0.5$n.ctrl.pdl) - 5
median(mc.res.0$n.ctrl.pdl) - 5
median(mc.res.0.5$n.ctrl.pdl) - 5



#####################
# Figure 2: performance of proposed estimator in terms of bias
# and efficiency
#####################
rsq.sum = rsq.mods[, .(bch.pdl.est = mean(bch.pdl.est),
                       pdl.est=mean(pdl.est),
                       offset.est=mean(pdl.offset.est),
                       allcontrol.est=mean(allcontrol.est),
                       bch.pdl.bias.perc=(mean(bch.pdl.est)-mean(alpha))/sd(bch.pdl.est),
                       pdl.bias.perc=(mean(pdl.est)-mean(alpha))/sd(pdl.est),
                       offset.bias.perc=(mean(pdl.offset.est)-mean(alpha))/sd(pdl.offset.est),
                       allcontrol.bias.perc=(mean(allcontrol.est)-mean(alpha))/sd(allcontrol.est)),
                   by=.(alpha.val=alpha, rsq)]

rsq.sum.long = melt(rsq.sum, 
                    id.vars = c("alpha.val", "rsq"), 
                    measure.vars = list("est"=3:6, "bias.perc"=7:10))
rsq.sum.long[,Method:=factor(variable,labels = c("Post-double Lasso (BCH penalty)", 
                                                 "Post-double Lasso (x-val penalty)", 
                                                 "Modified Post-double Lasso", 
                                                 "OLS with all controls"))]
rsq.sum.long[,alpha.label:=paste0('alpha*" = ',alpha.val,'"')]
rsq.bias.plot = ggplot(rsq.sum.long, aes(x=rsq, y=bias.perc, group=Method, shape=Method, 
                                         color=Method)) + 
  facet_wrap(~alpha.label, ncol = 2, labeller=label_parsed) + 
  geom_jitter(width=0.01) + 
  geom_hline(yintercept=0, linetype="dotted") + 
  geom_vline(xintercept=0) +
  theme_minimal() + 
  theme(legend.position="bottom", 
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"), 
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12)) + 
  scale_x_continuous(expand=c(0,0.02)) + 
  scale_color_manual(values=c("black", "black", "red", "black")) + 
  scale_y_continuous(breaks=c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5)) + 
  xlab("R squared") + 
  ylab("Bias as a fraction of empirical SD") + 
  guides(color = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2))

# efficiency analog
rsq.eff.sum = rsq.mods[, .(bch.pdl.sd = sd(bch.pdl.est),
                           pdl.sd=sd(pdl.est),
                           offset.sd=sd(pdl.offset.est),
                           allcontrol.sd=sd(allcontrol.est)),
                       by=.(alpha.val=alpha, rsq)]
rsq.eff.sum.long=melt(rsq.eff.sum, 
                      id.vars = c("alpha.val", "rsq"), 
                      measure.vars = list("sd"=3:6))
setnames(rsq.eff.sum.long, 4, "sd")
rsq.eff.sum.long[,Method:=factor(variable,labels = c("Post-double Lasso (BCH penalty)", 
                                                     "Post-double Lasso (x-val penalty)", 
                                                     "Modified Post-double Lasso", 
                                                     "OLS with all controls"))]
rsq.eff.sum.long[,alpha.label:=paste0('alpha*" = ',alpha.val,'"')]
rsq.eff.plot = ggplot(rsq.eff.sum.long, aes(x=rsq, y=sd, group=Method, shape=Method, 
                                            color=Method)) + 
  facet_wrap(~alpha.label, ncol = 2, labeller=label_parsed) + 
  geom_jitter(width=0.01) + 
  geom_hline(yintercept=0, linetype="dotted") + 
  geom_vline(xintercept=0) +
  theme_minimal() + 
  theme(legend.position="bottom", 
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"), 
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12)) + 
  scale_x_continuous(expand=c(0,0.02)) + 
  scale_color_manual(values=c("black", "black", "red", "black")) + 
  xlab("R squared") + 
  ylab("Empirical SD") + 
  guides(color = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2))

# Make Figure 2: Combined plot w/bias and efficiency
comb.plot = plot_grid(plot_grid(rsq.bias.plot + theme(legend.position = "none"), 
                                rsq.eff.plot + theme(legend.position = "none"),
                                
                                nrow = 2, 
                                ncol=1, 
                                labels = c("Bias","Standard Deviation"), 
                                label_x=0.54, label_y = 1, hjust = 0.5, scale=0.9),
                      get_legend(rsq.bias.plot + 
                                   theme(legend.position="right") + 
                                   guides(color = guide_legend(ncol = 1),
                                          shape = guide_legend(ncol = 1))), 
                      ncol=2, nrow=1, rel_widths = c(0.79, 0.21))
ggsave(here("figures","rsqcomb.png"), width=10, height=7)



#####################
# Figure 3: Summarize how alpha, p, and k, affect bias reduction 
#####################
alpha.dt = SummarizeMCRuns(alpha.mods)
np.dt = SummarizeMCRuns(p.mods)
k.dt = SummarizeMCRuns(k.mods)

# Does bias get appreciably worse or better as 
# the treatment effect is varied?
alpha.plot.data = alpha.dt[,.(alpha, 
                              bch.pdl.bias = (pdl.bch.est-alpha)/sd.bch.pdl.est,
                              xval.pdl.bias = (pdl.est-alpha)/sd.pdl.est,
                              offset.bias=(offset.est-alpha)/sd.offset.est),]

# As # confounders is varied
p.plot.data = np.dt[,.(n, p, 
                       bch.pdl.bias = (pdl.bch.est-alpha)/sd.bch.pdl.est,
                       xval.pdl.bias = (pdl.est-alpha)/sd.pdl.est,
                       offset.bias=(offset.est-alpha)/sd.offset.est),]

# As sparsity is varied
k.plot.data = k.dt[,.(k, 
                      bch.pdl.bias = (pdl.bch.est-alpha)/sd.bch.pdl.est,
                      xval.pdl.bias = (pdl.est-alpha)/sd.pdl.est,
                      offset.bias=(offset.est-alpha)/sd.offset.est),]


# Plot these
alpha.plot = ggplot(alpha.plot.data, aes(x=alpha)) + 
  geom_point(aes(y=bch.pdl.bias, 
                 shape="Post-double Lasso (BCH penalty)",
                 color="Post-double Lasso (BCH penalty)")) + 
  geom_point(aes(y=xval.pdl.bias, 
                 shape="Post-double Lasso (x-val penalty)",
                 color="Post-double Lasso (x-val penalty)")) + 
  geom_point(aes(y=offset.bias, 
                 shape="Modified post-double Lasso", 
                 color="Modified post-double Lasso", )) +
  geom_hline(yintercept=0, linetype="dotted") + 
  geom_vline(xintercept=0) +
  scale_color_manual(name="Method", values=c("red", "black", "black")) +  
  scale_shape_manual(name="Method", values=c(15,16,17)) +  
  scale_x_continuous(expand=c(0,0.1)) + 
  xlab(expression(alpha)) + 
  ylab("Bias as a fraction of Empirical SD") + 
  ggtitle(expression(alpha)) + 
  theme_minimal() + 
  theme(legend.position="bottom", 
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"), 
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5))

p.plot = ggplot(p.plot.data, aes(x=p)) + 
  geom_point(aes(y=bch.pdl.bias, 
                 shape="Post-double Lasso (BCH penalty)",
                 color="Post-double Lasso (BCH penalty)")) + 
  geom_point(aes(y=xval.pdl.bias, 
                 shape="Post-double Lasso (x-val penalty)",
                 color="Post-double Lasso (x-val penalty)")) + 
  geom_point(aes(y=offset.bias, 
                 shape="Modified post-double Lasso", 
                 color="Modified post-double Lasso", )) +
  geom_hline(yintercept=0, linetype="dotted") + 
  geom_vline(xintercept=0) +
  scale_color_manual(name="Method", values=c("red", "black", "black")) +  
  scale_shape_manual(name="Method", values=c(15,16,17)) +  
  scale_x_continuous(expand=c(0,2), 
                     breaks=seq(from=10, to=90, by=10),
                     labels = seq(from=10, to=90, by=10)) + 
  xlab("# candidate covariates (p)") + 
  ylab("Bias as a fraction of Empirical SD") + 
  ggtitle("p") + 
  theme_minimal() + 
  theme(legend.position="bottom", 
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"), 
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5))


k.plot = ggplot(k.plot.data, aes(x=k)) + 
  geom_point(aes(y=bch.pdl.bias, 
                 shape="Post-double Lasso (BCH penalty)",
                 color="Post-double Lasso (BCH penalty)")) + 
  geom_point(aes(y=xval.pdl.bias, 
                 shape="Post-double Lasso (x-val penalty)",
                 color="Post-double Lasso (x-val penalty)")) + 
  geom_point(aes(y=offset.bias, 
                 shape="Modified post-double Lasso", 
                 color="Modified post-double Lasso", )) +
  geom_hline(yintercept=0, linetype="dotted") + 
  geom_vline(xintercept=0) +
  scale_color_manual(name="Method", values=c("red", "black", "black")) +  
  scale_shape_manual(name="Method", values=c(15,16,17)) +  
  scale_x_continuous(expand = c(0,2), 
                     breaks = seq(from=0, to=25, by=5),
                     labels = seq(from=0, to=25, by=5)) + 
  xlab("Sparsity (k)") + 
  ylab("Bias as a fraction of Empirical SD") + 
  ggtitle("k") + 
  theme_minimal() + 
  theme(legend.position="bottom", 
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"), 
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5))

common.ylims = c(-0.6, 1.35)
noleg = theme(legend.position="none")
plot_grid(plot_grid(alpha.plot + ylim(common.ylims) + noleg,
                    p.plot + ylim(common.ylims) + noleg, 
                    k.plot + ylim(common.ylims) + noleg, 
                    ncol=3, align = T),
          get_legend(alpha.plot),
          ncol=1, nrow=2,
          rel_heights = c(0.9, 0.1)) +
  theme(plot.margin = unit(c(0.5,0.1,0.1,0.1), "cm"))
ggsave(here("figures","dgp_sens.png"), width=10, height=4)




#####################
# Make Table 1: inference
#####################
rsq.seq = unique(rsq.mods$rsq)
rsq.inf.sum = rsq.mods[, .(pdl.cover.auxadj = mean(pdl.offset.cover.auxadj),
                           pdl.cover.auxadjhc = mean(pdl.offset.cover.auxadjhc),
                           pdl.cover.noadj = mean(pdl.offset.cover),
                           pdl.offset.se.auxadj = mean(pdl.offset.se.auxadj),
                           pdl.offset.se.auxadjhc = mean(pdl.offset.se.auxadjhc),
                           pdl.offset.se.noadj = mean(pdl.offset.se),
                           pdl.offset.se.auxadj.fracols = mean(pdl.offset.se.auxadj)/mean(allcontrol.se),
                           pdl.offset.se.auxadjhc.fracols = mean(pdl.offset.se.auxadjhc)/mean(allcontrol.se),
                           pdl.offset.se.noadj.fracols = mean(pdl.offset.se)/mean(allcontrol.se),
                           pdl.offset.se.auxadj.fracolscjn = mean(pdl.offset.se.auxadj)/mean(allcontrol.se.cjn),
                           pdl.offset.se.auxadjhc.fracolscjn = mean(pdl.offset.se.auxadjhc)/mean(allcontrol.se.cjn),
                           pdl.offset.se.noadj.fracolscjn = mean(pdl.offset.se)/mean(allcontrol.se.cjn)
),
by=.(alpha.val=alpha, rsq)]

cover = dcast(data = rsq.inf.sum,
              formula = alpha.val ~ rsq, value.var = c("pdl.cover.auxadj"))[,`SE type`:="Homoskedastic"][,Stat:="Coverage"] 

cover.hc = dcast(data = rsq.inf.sum,
                 formula = alpha.val ~ rsq, value.var = c("pdl.cover.auxadjhc"))[,`SE type`:="HC"][,Stat:="Coverage"]

cover.noadj = dcast(data = rsq.inf.sum,
                    formula = alpha.val ~ rsq, value.var = c("pdl.cover.noadj"))[,`SE type`:="Unadjusted"][,Stat:="Coverage"]

se.frac = dcast(data = rsq.inf.sum,
                formula = alpha.val ~ rsq, value.var = c("pdl.offset.se.auxadj.fracolscjn"))[,`SE type`:="Homoskedastic"][,Stat:="SE (fraction of OLS SE)"] 

se.frac.hc = dcast(data = rsq.inf.sum,
                   formula = alpha.val ~ rsq, value.var = c("pdl.offset.se.auxadjhc.fracolscjn"))[,`SE type`:="HC"][,Stat:="SE (fraction of OLS SE)"]

se.frac.noadj = dcast(data = rsq.inf.sum,
                      formula = alpha.val ~ rsq, value.var = c("pdl.offset.se.noadj.fracolscjn"))[,`SE type`:="Unadjusted"][,Stat:="SE (fraction of OLS SE)"]


inf.table = rbindlist(list(cover,cover.hc, cover.noadj, se.frac, se.frac.hc, se.frac.noadj))
setcolorder(inf.table, c("alpha.val", "SE type", "Stat", as.character(rsq.seq)))
inf.table[,`SE type`:=factor(`SE type`, levels=c("Homoskedastic", "HC", "Unadjusted", ""), labels=c("Homoskedastic", "HC", "Unadjusted", ""))]
setorder(inf.table, alpha.val, `SE type`)
inf.table$`SE type`[c(2,4,6,8,10,12)] = "" # remove repeated labels for readability

inf.table.tex = kable(inf.table[,-1], format = "latex", digits = 3, booktabs = T,
                      caption="Coverage and standard error size for proposed adjusted standard errors. \\label{tab:se}")
inf.table.tex = group_rows(inf.table.tex, group_label="$\\mathbf{\\\\alpha=-1}$", start_row = 1, end_row = 6, escape = F, bold=F, latex_align="c")
inf.table.tex = group_rows(inf.table.tex, group_label="$\\mathbf{\\\\alpha=1}$", start_row = 7, end_row = 12, escape = F, bold=F, latex_align="c")
inf.table.tex = add_header_above(inf.table.tex, header=c(" "=2, "$R^2$"=7), escape = F)
inf.table.tex = add_footnote(inf.table.tex, 
                             label="Table notes: Coverage rates (nominal: 0.95) and standard errors as a fraction of those in Cattaneo, Jansson, \\& Newey (2018). Rows reflect the standard errors computed using (6) (``Homoskedastic'') and (6') (``HC''), and standard errors using OLS formulas (``Unadjusted''). Simulations are those underlying Figure 2.",
                             notation="none", escape=F, threeparttable = T)
kableExtra::save_kable(inf.table.tex, here("tables","setable.tex"))

