# Miscellaneous functions -------------------------------------------------------------------------------
#http://stackoverflow.com/questions/21011672/automatically-add-variable-names-to-elements-of-a-list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}

#function to fit and print a measurement model ------------------------------------------------------------
fitMM <- function (humspecmod,this1) {
  fithumspec <- sem(humspecmod, data=this1)
  summary(fithumspec, standardized=TRUE, fit.measures=TRUE)
  mod_ind <- modificationindices(fithumspec)
  print(head(mod_ind[order(mod_ind$mi, decreasing=TRUE), ], 4))
  return(fitMeasures(fithumspec))
}

#function to create and run a single indicator model and return omega and so on --------------------------
#  indicators are provided in character vector
#  pairwise correlations are optionally provided as a list of pairwise character vectors
singleindmodel <- function(indicators,paircorrs,dat) {
  pcs <- paste(unlist(lapply(paircorrs,FUN=function(x) paste(x[1], "~~", x[2], "\n"))),collapse=" ")
  SImod <- paste("factor =~", paste(indicators, collapse=" + "), "\n",
                 pcs,
                 "phantfac <~", paste(paste("1*",indicators,sep=""), collapse=" + "), "\n",
                 "factor ~~ 0*phantfac")

  fitSImod <- sem(SImod, dat[,indicators], estimator="ML")
  ParSImod <- parameterEstimates(fitSImod, standardized=TRUE)
  LoadingsSImod <- ParSImod[1:length(indicators), "std.all"]
  ErrorvarSImod <- 1 - LoadingsSImod^2
  ImpliedCorrSImod <- lavTech(fitSImod, what='cor.lv')
  OmegaSImod <- ImpliedCorrSImod[[1]][1,2]^2            #Squared correlation between factor and phantom variable

  varSImod <- var(apply(dat[,indicators], MARGIN=1, FUN=mean)) #variance of composite
  SDSImod <- sqrt(varSImod)               #SD of composite, to pass back for analysis
  eSImod <- (1-OmegaSImod)*varSImod       #error term of single indicator
  return(listN(OmegaSImod,eSImod,SDSImod))
}
