# Code to run Monte Carlo simulations underlying main
# results. 
# Saves simulation results to .Rdata files in results directory.
# The separate SummarizeResults.R file takes those .Rdata files
# and produces figures, a table, and key summary statistics.
#
# Note: in each group of simulations (baseline, varying R^2, varying alpha,
# varying p, and varying k), the random seed is set separately
# (to a different fixed value). This is not strictly necessary, but
# is done to facilitate running a specific simulation group 
# in isolation while maintaining reproducibility.
#
# Author: Steve Miller (steve.j.miller@colorado.edu)

# Load libraries
library(here)
library(MASS)
library(fixest)
library(data.table)
library(glmnet)
library(hdm)
library(sandwich)

# Global: # Monte Carlo replicates
n.mc.rep = 1000

# Functions for SEs in high dimensional regression 
# adapted from Cattaneo, Jansson, & Newey (2018) replication files:
# https://github.com/mdcattaneo/replication-CJN_2018_JASA/blob/master/main.fun.R
XXinv = function(x) chol2inv(chol(crossprod(x)))
CJNSE = function(Y,X,W) {
  # Compute M, kappa
  K = ncol(W); 
  qr.W = qr(W); 
  if (qr.W$rank==K){
    M = -tcrossprod(qr.Q(qr.W))
    diag(M) = 1 + diag(M)
    kappa = try(chol2inv(chol(M*M)), TRUE);
  } else {
    error("W not full rank")
  }
  # Form components of SE estimate
  MX = M%*%X 
  XMX_inv = XXinv(MX) 
  XMY = crossprod(MX,Y) 
  n = length(Y)
  beta.hat = XMX_inv%*%XMY
  u2.hat = (M%*%Y - MX%*%beta.hat)^2
  
  # Form variance estimate
  if (is.matrix(kappa)){
    VHK = XMX_inv * crossprod(MX*(kappa%*%u2.hat),MX) * XMX_inv
  } else {VHK = NA}
  
  # Give back SE
  return(sqrt(diag(VHK)))
}
  

# Main function to run a Monte Carlo simulation
RunMC = function(n=500, 
                 k=5, 
                 p=200,
                 n.replicates=1000,
                 alpha=0,
                 confounding.strength=0.5,
                 sigma.y=1,
                 sigma.d=1,
                 run.bch.pdl = T,
                 compute.cjnse = T,
                 vcov.offset = "iid",
                 trace=F) {
  
  beta = c(rep(confounding.strength,k), rep(0,p-k))
  gamma = c(rep(confounding.strength,k), rep(0,p-k))
  
  all.est = rbindlist(lapply(1:n.replicates, FUN=function(replicate) {
    if(trace) {
      print(paste0("replicate ", replicate, " of ", n.replicates))
    }
    X = mvrnorm(n = n, mu = rep(0, p), Sigma = diag(p))
    nu = rnorm(n = n, mean = 0, sd = 1)
    eta = rnorm(n = n, mean = 0, sd = 1)
    D = X %*% gamma + sigma.d*nu
    if(is.numeric(sigma.y)) {
      Y = D*alpha + X %*% beta + sigma.y*eta
    } else {
      if(sigma.y == "X1") {
        Y = D*alpha + X %*% beta + abs(X[,1])*eta
      }
    }
    replicate.data = data.table(cbind(Y, D, X))
    setnames(replicate.data, c("Y", "D", paste0("X",1:p)))
    

    mod.allcontrol = feols(as.formula(paste0("Y ~ 0 + D + ", 
                                             paste0("X",1:p, collapse= " + "))), data = replicate.data)
    if(compute.cjnse) {
      allcontrol.se.cjn = CJNSE(Y=replicate.data$Y,
                                X=as.matrix(x=replicate.data$D,ncol=1),
                                W=X
      )
    } else {
      allcontrol.se.cjn = NA
    }
    
    
    # Find controls from selection/treatment equation
    pdl.selection.cvres = cv.glmnet(x=X,
                                    y=D,
                                    family="gaussian",
                                    alpha=1
                                    )
    selection.controls = which(as.numeric(coef(pdl.selection.cvres, s="lambda.min"))[-1] != 0)
    n.corr.sel.controls = sum(selection.controls %in% 1:k)

    # Find controls from outcome equation
    pdl.outcome.cvres = cv.glmnet(x=X,
                                 y=Y,
                                 family="gaussian",
                                 alpha=1)
    outcome.controls = which(as.numeric(coef(pdl.outcome.cvres, s="lambda.min"))[-1] != 0)

    # Combine controls and estimate
    all.controls = union(selection.controls, outcome.controls)
    n.corr.pdl.controls = sum(all.controls %in% 1:k)
    if(length(all.controls)>0) {
      mod.pdl = feols(as.formula(paste0("Y ~ 0 + D + ", 
                                 paste0("X",all.controls, collapse= " + "))), data = replicate.data)
    } else {
      mod.pdl = feols(Y ~ 0 + D, data = replicate.data)
    }
    
    # Post-double Lasso as described in Belloni, Chernozhukov, & Hansen (2014)
    # with theoretically-motivated penalty parameter
    if(run.bch.pdl) {
      bch.pdl = hdm::rlassoEffect(x=X, y=Y, d=D)
      bch.pdl.est = bch.pdl$coefficients
      bch.pdl.se = bch.pdl$se
      bch.pdl.cover = alpha >= (bch.pdl.est - 1.96*bch.pdl.se) & alpha <= (bch.pdl.est + 1.96*bch.pdl.se)
    } else {
      bch.pdl.est = NA
      bch.pdl.se = NA
      bch.pdl.cover = NA
    }
    
    X.withtreatment = cbind(D, X)
    
    # Proposed Estimator -- get OLS estimate of treatment effect, include that
    # as an offset in lasso for outcome equation (equivalent to removing from Y) 

    pre.ols.offset = mod.allcontrol$coefficients["D"]*D
    pdl.outcome.offset.cvres = cv.glmnet(x=X,
                                                y=Y,
                                                offset = pre.ols.offset,
                                                family="gaussian",
                                                alpha=1)
    outcome.offset.controls = which(as.numeric(coef(pdl.outcome.offset.cvres, s="lambda.min"))[-1] != 0)
    
    # Combine controls & estimate
    all.offset.controls = union(selection.controls, outcome.offset.controls)
    n.corr.offset.controls = sum(all.offset.controls %in% 1:k)
    
    if(length(all.offset.controls)>0) {
      if(vcov.offset=="iid") {
        mod.offset.pdl = feols(as.formula(paste0("Y ~ 0 + D + ", 
                                                      paste0("X",all.offset.controls, collapse= " + "))), 
                               vcov=vcov.offset,
                               data = replicate.data)
      } else if(vcov.offset %in% c("HC0","HC1","HC2","HC3","HC4")) {
        mod.offset.pdl = feols(as.formula(paste0("Y ~ 0 + D + ", 
                                                 paste0("X",all.offset.controls, collapse= " + "))), 
                               vcov = function(x) sandwich::vcovHC(x, type = vcov.offset),
                               data = replicate.data)
      }
      if(compute.cjnse) {
      # Compute Cattaneo, Jansson, & Newey (2018) SEs for proposed model
        pdl.offset.se.cjn = CJNSE(Y=replicate.data$Y,
                                  X=as.matrix(x=replicate.data$D,ncol=1),
                                  W=as.matrix(replicate.data[,paste0("X",all.offset.controls),with=F],
                                              ncol=length(all.offset.controls))
        )
      } else {
        pdl.offset.se.cjn = NA
      }
    } else {
      if(vcov.offset=="iid") {
        mod.offset.pdl = feols(Y ~ 0 + D, 
                               vcov=vcov.offset,
                               data = replicate.data)
      } else if(vcov.offset %in% c("HC0","HC1","HC2","HC3","HC4")) {
        mod.offset.pdl = feols(Y ~ 0 + D, 
                               vcov = function(x) sandwich::vcovHC(x, type = vcov.offset),
                               data = replicate.data)
      }
      pdl.offset.se.cjn = NA
    }
    
    # The estimate of sigma^2 with this approach is likely to be downward biased (Reid et al., 2016)
    # resulting in SEs that are too small.
    # One alternative is to estimate sigma^2 by an auxiliary cross-validated 
    # Lasso fit of the structural model
    aux.lasso.cvres = cv.glmnet(x=X.withtreatment,
                                y=Y,
                                family="gaussian",
                                alpha=1) 
    # calculate s_{hat}, including intercept
    n.sel.coef = length(which(as.numeric(coef(aux.lasso.cvres, s="lambda.min")) !=0))
    aux.lasso.sq.resid = (as.numeric(Y-predict(aux.lasso.cvres, newx = X.withtreatment, s = "lambda.min")))^2
    sigma.sq.aux.lasso = 1/(n-n.sel.coef)*sum(aux.lasso.sq.resid)
    pdl.offset.se.auxadj = sigma.sq.aux.lasso/mod.offset.pdl$sigma2*coeftable(mod.offset.pdl)[1,2]
    
    # HC3 variant using lasso residuals
    if(length(all.offset.controls)> 0) {
      X.offset = as.matrix(replicate.data[,c("D",paste0("X",all.offset.controls)),with=F])  
    } else {
      X.offset = as.matrix(replicate.data[,"D",with=F])  
    }
    X.offset.prime = t(X.offset)
    bread=XXinv(X.offset)
    hatmat = X.offset %*% bread %*% X.offset.prime
    meat=X.offset.prime %*% diag(x=aux.lasso.sq.resid/(1-diag(hatmat)), nrow = n, ncol=n) %*% X.offset
    pdl.offset.se.auxadjhc = sqrt(diag(bread %*% meat %*% bread))[1]

    # assemble results
    estimates = data.table(# Document some of the parameters supplied
                           alpha = alpha,
                           n=n, 
                           k=k, 
                           p=p,
                           confounding.strength=confounding.strength,
                           sigma.y=sigma.y,
                           sigma.d=sigma.d,
                           # OLS on full set of controls
                           allcontrol.est=coef(mod.allcontrol)["D"],
                           allcontrol.se=coeftable(mod.allcontrol)[1,2],
                           allcontrol.se.cjn = allcontrol.se.cjn,
                           allcontrol.cover=alpha >= (coef(mod.allcontrol)["D"] - 1.96*coeftable(mod.allcontrol)[1,2]) & alpha <= (coef(mod.allcontrol)["D"] + 1.96*coeftable(mod.allcontrol)[1,2]),
                           allcontrol.cover.cjn = alpha >= (coef(mod.allcontrol)["D"] - 1.96*allcontrol.se.cjn) & alpha <= (coef(mod.allcontrol)["D"] + 1.96*allcontrol.se.cjn),

                           # Point estimates for various PDL estimators
                           bch.pdl.est = bch.pdl.est,
                           pdl.est=coef(mod.pdl)["D"],
                           pdl.offset.est = coef(mod.offset.pdl)["D"],
                           
                           # SEs for PDL estimators
                           bch.pdl.se = bch.pdl.se,
                           pdl.se = coeftable(mod.pdl)[1,2],
                           pdl.offset.se = coeftable(mod.offset.pdl)[1,2],
                           pdl.offset.se.auxadj = pdl.offset.se.auxadj,
                           pdl.offset.se.auxadjhc = pdl.offset.se.auxadjhc,
                           pdl.offset.se.cjn = pdl.offset.se.cjn,

                           # Coverage for PDL estimators
                           bch.pdl.cover = bch.pdl.cover,
                           pdl.cover = alpha >= (coef(mod.pdl)["D"] - 1.96*coeftable(mod.pdl)[1,2]) & alpha <= (coef(mod.pdl)["D"] + 1.96*coeftable(mod.pdl)[1,2]),
                           pdl.offset.cover = alpha >= (coef(mod.offset.pdl)["D"] - 1.96*coeftable(mod.offset.pdl)[1,2]) & alpha <= (coef(mod.offset.pdl)["D"] + 1.96*coeftable(mod.offset.pdl)[1,2]),
                           pdl.offset.cover.auxadj = alpha >= (coef(mod.offset.pdl)["D"] - 1.96*pdl.offset.se.auxadj) & alpha <= (coef(mod.offset.pdl)["D"] + 1.96*pdl.offset.se.auxadj),
                           pdl.offset.cover.auxadjhc = alpha >= (coef(mod.offset.pdl)["D"] - 1.96*pdl.offset.se.auxadjhc) & alpha <= (coef(mod.offset.pdl)["D"] + 1.96*pdl.offset.se.auxadjhc),
                           pdl.offset.cover.cjn = alpha >= (coef(mod.offset.pdl)["D"] - 1.96*pdl.offset.se.cjn) & alpha <= (coef(mod.offset.pdl)["D"] + 1.96*pdl.offset.se.cjn),

                           # Some selection diagnostics for cross-validated 
                           # PDL estimators
                           n.ctrl.pdl = length(all.controls),
                           n.ctrl.offset = length(all.offset.controls),
                           n.corr.sel.controls = n.corr.sel.controls,
                           n.corr.pdl.controls = n.corr.pdl.controls,
                           n.corr.offset.controls = n.corr.offset.controls
    )
    return(estimates)
  }))
}


##############################
# Run simulations
##############################


# Baseline MC simulations to illustrate bias
set.seed(0)
mc.res._0.5 = RunMC(alpha=-0.5, trace=T)
mc.res.0 = RunMC(alpha=0, trace=T)
mc.res.0.5 = RunMC(alpha=0.5, trace=T)
save(mc.res._0.5, mc.res.0, mc.res.0.5, file=here("results","simpleres.Rdata"))

# Assess how bias of existing and proposed method perform 
# as DGP characteristics are varied.

# Confounding strength. Use the same setup as those underlying 
# Wuthrich & Zhu (2023) DGP A5 & A6, where regular post double Lasso struggles
# even with cross-validated lambda
set.seed(pi)
rsq.seq = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
n.obs = 500
p = 200
k = 5
var.x = 1
alphas = c(-1, 1)

rsq.mods = rbindlist(lapply(alphas, FUN=function(alpha) {
  rbindlist(lapply(rsq.seq, FUN=function(rsq) {
    confounding.strength = sqrt(rsq/(var.x*k*(1 - rsq)))
    print(paste(alpha, rsq, confounding.strength))
    res = RunMC(alpha=alpha, 
                n=n.obs, 
                n.replicates=n.mc.rep,
                confounding.strength = confounding.strength, 
                p = p, 
                trace=T)
    res$rsq = rsq
    return(res)
  }))
}))

save(rsq.mods, file=here("results","rsqres.Rdata"))

#########################################
# Magnitude of treatment effect (alpha)
#########################################
set.seed(42)
alpha.seq = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5)
rsq = 0.1
k = 5
var.x = 1
confounding.strength = sqrt(rsq/(var.x*k*(1 - rsq)))
n.obs = 200
p = 100

alpha.mods = lapply(alpha.seq, FUN=function(alpha) {
  print(alpha)
  res = RunMC(alpha=alpha, n=n.obs, n.replicates=n.mc.rep, p = p, confounding.strength = confounding.strength, trace=T)
  return(res)
})
save(alpha.mods, file=here("results","alphares.Rdata"))

#########################################
# Number of candidate confounders (p)
#########################################
set.seed(8675309)
alpha = 1
n.obs = 200
rsq = 0.1
k = 5
var.x = 1
confounding.strength = sqrt(rsq/(var.x*k*(1 - rsq)))
p.seq = seq(from=10, to=90, by=10)
p.mods = lapply(p.seq, FUN=function(p) {
  print(p)
  p.res = RunMC(alpha=alpha, n=n.obs, n.replicates=n.mc.rep, p = p, trace=T)
  return(p.res)
})
save(p.mods, file=here("results","pres.Rdata"))

#########################################
# Sparsity (k)
#########################################
set.seed(11)
alpha = 1
n.obs = 200
p = 100
rsq = 0.1
var.x = 1
k.seq = c(1, 5, 10, 15, 20, 25)

k.mods = lapply(k.seq, FUN=function(k) {
  print(k)
  confounding.strength = sqrt(rsq/(var.x*k*(1 - rsq))) # vary confounding strength to hold R^2 fixed
  k.res = RunMC(alpha=alpha, n=n.obs, n.replicates=n.mc.rep, p = p, k = k, confounding.strength = confounding.strength, trace=T)
  return(k.res)
})
save(k.mods, file=here("results","kres.Rdata"))

