library(tidyr)
library(TMB)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)

# Setup Data ####
source("src/load_format_data.R")
source("src/Extract_model_fits.R")
sd(temp.data$temp_na_rm)
sd(seal.data$count)

temp_obs_prior <- c(log(0.25), 0.05)
seal_obs_prior <- c(log(0.25), 0.05)

# state space contracaecum only - used for initialization ####
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

base_ests <- extract.fit(rep, yearlist = yearlist)

setup_tmb_data <- function(base.ests,
                           data1,
                           data2,
                           yearindex1,
                           yearindex2,
                           yearlist,
                           ndata,
                           Xij,
                           Uij,
                           y,
                           data1_obs_prior,
                           data2_obs_prior) {
  
  x1tobs <- data1
  x1yearindex <- yearindex1
  colnames(x1tobs) <- c("year", "xt")
  x1tobs$xt <- as.numeric(scale(x1tobs$xt))
  nx1 <- nrow(x1tobs)
  x2tobs <- data2
  x2yearindex <- yearindex2
  colnames(x2tobs) <- c("year", "xt")
  x2tobs$xt <- as.numeric(scale(x2tobs$xt))
  nx2 <- nrow(x2tobs)
  
  
  tmb_data <- list(
    n = ndata,
    nyears = length(yearlist),
    y = y,
    Xij = Xij,
    Uij = Uij,
    yearindex = yearindex - 1,
    x1tobs = x1tobs$xt,
    x1yearindex = x1yearindex - 1,
    nx1 = nx1,
    x2tobs = x2tobs$xt,
    x2yearindex = x2yearindex - 1,
    nx2 = nx2,
    priorlogsigmaobs1_mu = data1_obs_prior[1],
    priorlogsigmaobs1_sigma = data1_obs_prior[2],
    priorlogsigmaobs2_mu = data2_obs_prior[1],
    priorlogsigmaobs2_sigma = data2_obs_prior[2]
  )
  return(tmb_data)
}


# Full Model ####
tmb_data <- setup_tmb_data(base.ests = base_ests,
                           data1 = seal.data,
                           data2 = temp.data,
                           yearindex1 = syearindex,
                           yearindex2 = tyearindex,
                           yearlist = yearlist,
                           ndata = ndata,
                           Xij = Xij,
                           Uij = Uij,
                           y = spc.data$count,
                           data1_obs_prior = seal_obs_prior,
                           data2_obs_prior = temp_obs_prior)
x1obs[,2] <- tmb_data$x1tobs
x2obs[,2] <- tmb_data$x2tobs
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
                  logsigma_x1obs = seal_obs_prior[1],
                  x1o = 0,
                  thetax1 = 0,
                  vx2t = rep(0, length(yearlist) - 1),
                  logsigma_vx2t = 0,
                  logsigma_x2obs = temp_obs_prior[1],
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
    random = c("gamma", "vt", "vx1t", "vx2t", "logsigma_x1obs", "logsigma_x2obs"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)

both_ests <- extract.fit(rep, yearlist = yearlist, plotx2t = T, x1tobs = x1obs, x2tobs = x2obs)

model.output <- list(obj = obj, opt = opt, rep = rep)
save(file = "src/both_predictors.rda", model.output)

# Seal Model ####
tmb_data <- setup_tmb_data(base.ests = base_ests,
                           data1 = seal.data,
                           data2 = temp.data,
                           yearindex1 = syearindex,
                           yearindex2 = tyearindex,
                           yearlist = yearlist,
                           ndata = ndata,
                           Xij = Xij,
                           Uij = Uij,
                           y = spc.data$count,
                           data1_obs_prior = seal_obs_prior,
                           data2_obs_prior = temp_obs_prior)
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
                  logsigma_x1obs = seal_obs_prior[1],
                  x1o = 0,
                  thetax1 = 0,
                  vx2t = rep(0, length(yearlist) - 1),
                  logsigma_vx2t = 0,
                  logsigma_x2obs = temp_obs_prior[1],
                  x2o = 0)



model <- "one_covariate_state_space"
compile(paste0("src/TMB/", model, ".cpp"))

dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt", "vx1t", "vx2t",  "logsigma_x1obs", "logsigma_x2obs"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
seal_ests <- extract.fit(rep, yearlist = yearlist, plotx2t = T, x1tobs = x1obs, x2tobs = x2obs)

model.output <- list(obj = obj, opt = opt, rep = rep)
save(file = "src/seal_predictors.rda", model.output)

# Temp Model ####
tmb_data <- setup_tmb_data(base.ests = base_ests,
                           data1 = temp.data,
                           data2 = seal.data,
                           yearindex1 = tyearindex,
                           yearindex2 = syearindex,
                           yearlist = yearlist,
                           ndata = ndata,
                           Xij = Xij,
                           Uij = Uij,
                           y = spc.data$count,
                           data1_obs_prior = temp_obs_prior,
                           data2_obs_prior = seal_obs_prior)
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
                  logsigma_x1obs = temp_obs_prior[1],
                  x1o = 0,
                  thetax1 = 0,
                  vx2t = rep(0, length(yearlist) - 1),
                  logsigma_vx2t = 0,
                  logsigma_x2obs = seal_obs_prior[1],
                  x2o = 0)



model <- "one_covariate_state_space"
compile(paste0("src/TMB/", model, ".cpp"))

dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt", "vx1t", "vx2t",  "logsigma_x1obs", "logsigma_x2obs"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
temp.ests <- extract.fit(rep, x1tobs = x2obs, x2tobs = x1obs, yearlist = yearlist, plotx1t = T, plotx2t = T)

model.output <- list(obj = obj, opt = opt, rep = rep)
save(file = "src/temp_predictors.rda", model.output)

# No predictors ####
tmb_data <- setup_tmb_data(base.ests = base_ests,
                           data1 = seal.data,
                           data2 = temp.data,
                           yearindex1 = syearindex,
                           yearindex2 = tyearindex,
                           yearlist = yearlist,
                           ndata = ndata,
                           Xij = Xij,
                           Uij = Uij,
                           y = spc.data$count,
                           data1_obs_prior = seal_obs_prior,
                           data2_obs_prior = temp_obs_prior)
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
                  logsigma_x1obs = seal_obs_prior[1],
                  x1o = 0,
                  vx2t = rep(0, length(yearlist) - 1),
                  logsigma_vx2t = 0,
                  logsigma_x2obs = temp_obs_prior[1],
                  x2o = 0)


model <- "no_covariate_state_space"
compile(paste0("src/TMB/", model, ".cpp"))

dyn.load(dynlib(paste0("src/TMB/",model)))
obj <-
  MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    DLL = model,
    random = c("gamma", "vt", "vx1t", "vx2t",  "logsigma_x1obs", "logsigma_x2obs"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
extract.fit(rep, x1tobs = x1obs, x2tobs = x2obs, yearlist = yearlist, plotx1t = T, plotx2t = T)
model.output <- list(obj = obj, opt = opt, rep = rep)
save(file = "src/no_predictors.rda", model.output)

# Model Selection ####
TMBAIC=function(opt, p=2, n=Inf){
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
  return( Return )
}

filenames <- c("both_predictors", "seal_predictors", "temp_predictors", "no_predictors")

AIC <- matrix(nrow = 4, ncol = 2)
rownames(AIC) <- filenames

for (i in 1:length(filenames)) {
  filename <- paste0("src/", filenames[i], ".rda")
  load(filename)
  AIC[i,1] <- TMBAIC(model.output$opt, n = ndata)
  
}
AIC[,2] <- AIC[,1] - min(AIC[,1])
colnames(AIC) <- c("AIC", "delta AIC")
print(AIC)
