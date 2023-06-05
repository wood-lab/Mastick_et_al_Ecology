extract.fit <-
  function(rep,
           yearlist,
           x1tobs = NULL,
           x2tobs = NULL,
           plotnt = T,
           plotrt = T,
           plotx1t = F,
           plotx2t = F) {
    # extract stuff and make plots
    fixef <- summary(rep, "fixed")
    raef <- summary(rep, "random")
    reef <- summary(rep, "report")
    
    get_est <- function(ef, parname)    return(ef[grep(rownames(ef), pattern = paste0("\\b", parname, "\\b")), ])
    
    beta_est <- get_est(fixef, "beta")
    #### Extract random effects ####
    gamma_est <- get_est(raef, "gamma")
    
    #### Extract population size ####
    # plot nt
    if (plotnt) {
      ntEst <- get_est(reef, "lognt")
      
      plotdf <- tibble(
        year = yearlist,
        Nt = exp(ntEst[, 1]),
        Ntlower = exp(ntEst[, 1] - ntEst[, 2]),
        Ntupper = exp(ntEst[, 1] + ntEst[, 2])
      )
      ntplot <- ggplot(data = plotdf, aes(x = year, y = Nt)) +
        geom_ribbon(aes(x = year, ymin = Ntlower, ymax = Ntupper), fill = "gray") +
        geom_line(linewidth=1.5)+theme_classic()+ theme(text = element_text(size = 14))+
        xlab("Year")+ylab("Estimated Contracaecum")+
        scale_x_continuous(limits=c(1900, 2020))+ggtitle('E')
      print(ntplot)
    }
    
    if (plotrt) {
      rtEst <- get_est(reef, "rt")
      plotdf <- tibble(
        year = yearlist[-length(yearlist)],
        rt = rtEst[,1],
        rtlower = rtEst[, 1] - rtEst[, 2],
        rtupper = rtEst[, 1] + rtEst[, 2]
      )
      rtplot <- ggplot(data = plotdf, aes(x = year, y = rt)) +
        geom_ribbon(aes(x = year, ymin = rtlower, ymax = rtupper), fill = "gray") +
        geom_line()+theme_classic()+ theme(text = element_text(size = 14))+xlab("Year")
      print(rtplot)
    }
    
    plotxtfun <- function(reef, var, yearlist, xtobs) {
      xtEst <- get_est(reef, var)
      
      plotdf <- tibble(
        year = yearlist,
        Xt = xtEst[, 1],
        Xtlower = xtEst[, 1] - xtEst[, 2],
        Xtupper = xtEst[, 1] + xtEst[, 2]
      )
      xtplot <- ggplot(data = plotdf, aes(x = year, y = Xt)) +
        geom_ribbon(aes(x = year, ymin = Xtlower, ymax = Xtupper), fill = "gray") +
        geom_line() +
        geom_point(data = xtobs, aes(x = year, y = xt))+theme_classic()+ theme(text = element_text(size = 14))+xlab("Year")
      print(xtplot)
      
      
    }
    
    if (plotx1t) {
      plotxtfun(reef, "x1t", yearlist, xtobs = x1tobs)
      thetax <- get_est(fixef, "thetax")
      print(thetax)
    }
    
    if (plotx2t) {
      plotxtfun(reef, "x1t", yearlist, xtobs = x1tobs)
      thetax1 <- get_est(fixef, "thetax1")
      print(thetax1)
      plotxtfun(reef, "x2t", yearlist, xtobs = x2tobs)
      thetax2 <- get_est(fixef, "thetax2")
      print(thetax2)
      
    }
    # get vt's for model fitting
    vtest <- get_est(raef, "vt")
    logsigma_vt <- get_est(fixef, "logsigma_vt")
    logitrho <- get_est(fixef, "logitrho")
    logsigma_gamma <- get_est(fixef, "logsigma_gamma")
    logno <- get_est(fixef, "logno")
    logphi <- get_est(fixef, "logphi")
    
    return(list (beta = as.numeric(beta_est[,1]),
                 gamma = as.numeric(gamma_est[,1]),
                 vt = as.numeric(vtest[,1]),
                 logsigma_vt = as.numeric(logsigma_vt[1]),
                 logitrho = as.numeric(logitrho[1]),
                 logsigma_gamma = as.numeric(logsigma_gamma[1]),
                 logno = as.numeric(logno[1]),
                 logphi = as.numeric(logphi[1])
    )
    )
    
  }
