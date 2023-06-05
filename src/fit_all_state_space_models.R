library(tidyr)
library(TMB)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)

# Setup Data ####
source("src/load_format_data.R")
source("src/Extract_model_fits.R")

pdf(file = "src/summary_plots.pdf",
    height =8,
    width = 8)
# no Covariates ####
tmb_data <- list(n = ndata,
                 nyears = nyears,
                 y = spc.data$count,
                 Xij = Xij,
                 Uij = Uij,
                 yearindex = yearindex - 1
)

tmb_pars <- list( beta = rep(0, 3),
                  gamma = rep(0, length(host.2.keep)),
                  logsigma_gamma = 0,
                  vt = rep(0, nyears -1),
                  logitrho = 0,
                  logsigma_vt = 0,
                  logno = 0,
                  logphi = 0
)
model <- "state_space"
compile(paste0("src/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
model.output <- list(obj = obj, opt = opt, rep = rep)
base_ests <- extract.fit(rep, yearlist = yearlist)

save(file = "src/no_covars.rda", model.output)

# Seal Only ####
xtobs <- seal.data
xyearindex <- syearindex 
colnames(xtobs) <- c("year", "xt")
xtobs$xt <- as.numeric(scale(xtobs$xt))
nx <- nrow(xtobs)

tmb_data <- list(n = ndata,
                 nyears = length(yearlist),
                 y = spc.data$count,
                 Xij = Xij,
                 Uij = Uij,
                 yearindex = yearindex - 1,
                 xtobs = xtobs$xt,
                 xyearindex = xyearindex - 1,
                 nx = nx
)


tmb_pars <- list( beta = base_ests$beta,
                  gamma = base_ests$gamma,
                  logsigma_gamma = base_ests$logsigma_gamma,
                  vt = base_ests$vt,
                  logitrho = base_ests$logitrho,
                  logsigma_vt = base_ests$logsigma_vt,
                  logno = base_ests$logno,
                  logphi = base_ests$logphi,
                  vxt = rep(0, length(yearlist) - 1),
                  logsigma_vxt = 0,
                  logsigma_xobs = -.5,
                  xo = 0,
                  thetax = 0)

model <- "covariate_state_space"
compile(paste0("src/TMB/", model, ".cpp"))

dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt", "vxt"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)

seal_ests <- extract.fit(rep, yearlist = yearlist, plotxt = T, xtobs = xtobs)

model.output <- list(obj = obj, opt = opt, rep = rep)
save(file = "src/seals.rda", model.output)

# Temperature Only ####
xtobs <- temp.data
xyearindex <- tyearindex 
colnames(xtobs) <- c("year", "xt")
xtobs$xt <- as.numeric(scale(xtobs$xt))
nx <- nrow(xtobs)

tmb_data <- list(n = ndata,
                 nyears = length(yearlist),
                 y = spc.data$count,
                 Xij = Xij,
                 Uij = Uij,
                 yearindex = yearindex - 1,
                 xtobs = xtobs$xt,
                 xyearindex = xyearindex - 1,
                 nx = nx,
                 priorlogsigmaobs_mu = log(0.25),
                 priorlogsigmaobs_sigma = 0.25
)


tmb_pars <- list( beta = base_ests$beta,
                  gamma = base_ests$gamma,
                  logsigma_gamma = base_ests$logsigma_gamma,
                  vt = base_ests$vt,
                  logitrho = base_ests$logitrho,
                  logsigma_vt = base_ests$logsigma_vt,
                  logno = base_ests$logno,
                  logphi = base_ests$logphi,
                  vxt = rep(0, length(yearlist) - 1),
                  logsigma_vxt = 0,
                  logsigma_xobs = -.5,
                  xo = 0,
                  thetax = 0)

model <- "covariate_state_space"
compile(paste0("src/TMB/", model, ".cpp"))

dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt", "vxt"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)

temp_ests <- extract.fit(rep, yearlist = yearlist, plotxt = T, xtobs = xtobs)

model.output <- list(obj = obj, opt = opt, rep = rep)
save(file = "src/temp.rda", model.output)

# Full Model ####
x1tobs <- temp.data
x1yearindex <- tyearindex 
colnames(x1tobs) <- c("year", "xt")
x1tobs$xt <- as.numeric(scale(x1tobs$xt))
nx1 <- nrow(x1tobs)
x2tobs <- seal.data
x2yearindex <- syearindex 
colnames(x2tobs) <- c("year", "xt")
x2tobs$xt <- as.numeric(scale(x2tobs$xt))
nx2 <- nrow(x2tobs)


tmb_data <- list(n = ndata,
                 nyears = length(yearlist),
                 y = spc.data$count,
                 Xij = Xij,
                 Uij = Uij,
                 yearindex = yearindex - 1,
                 x1tobs = x1tobs$xt,
                 x1yearindex = x1yearindex - 1,
                 nx1 = nx1,
                 x2tobs = x2tobs$xt,
                 x2yearindex = x2yearindex - 1,
                 nx2 = nx2
)


tmb_pars <- list( beta = base_ests$beta,
                  gamma = base_ests$gamma,
                  logsigma_gamma = base_ests$logsigma_gamma,
                  vt = base_ests$vt,
                  logitrho = base_ests$logitrho,
                  logsigma_vt = base_ests$logsigma_vt,
                  logno = base_ests$logno,
                  logphi = base_ests$logphi,
                  vx1t = rep(0, length(yearlist) - 1),
                  logsigma_vx1t = 0,
                  logsigma_x1obs = -.5,
                  x1o = 0,
                  thetax1 = 0,
                  vx2t = rep(0, length(yearlist) - 1),
                  logsigma_vx2t = 0,
                  logsigma_x2obs = -.5,
                  x2o = 0,
                  thetax2 = 0)

model <- "two_covariate_state_space"
compile(paste0("src/TMB/", model, ".cpp"))

dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt", "vx1t", "vx2t"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)

both_ests <- extract.fit(rep, yearlist = yearlist, plotx2t = T, xtobs = x1tobs, x2tobs = x2tobs)

model.output <- list(obj = obj, opt = opt, rep = rep)
save(file = "src/both.rda", model.output)
dev.off()
