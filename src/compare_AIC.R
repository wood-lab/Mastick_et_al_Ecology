# Compare Models
model.names <- c("no_covars", "temp", "seals", "temp_and_seals")

TMBAIC=function(opt, p=2, n=Inf){
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
  return( Return )
}


AIC <- matrix(NA, nrow= length(model.names), ncol = 2)
rownames(AIC) <- model.names

for (i in 1:length(model.names)) {
  filename <- paste0("src/",model.names[i],".rda")
  load(file = filename)
  AIC[i,1] <- TMBAIC(model.output$opt)
}
