library(here)
library(parallel)
library(coda)
library(rcarbon)
library(nimbleCarbon)
load(here('c14dates','02_rice_data.RData'))

#MCMC settings ----
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12,34,56,78)
niter  <- 600000
nburnin  <- 300000
thin  <- 30

# Run Model
out  <-  parLapply(cl = cl, X = seeds, fun = runscript, d = dat.rice, constants = constants.rice, inits = inits.rice, niter = niter, nburnin = nburnin,thin = thin)
stop(cl)

rhats.rice <- gelman.diag(out)
rhats.rice[[1]][which(rhats.rice[[1]][,1]>1.01),]
ess.rice <- effectiveSize(out)
post.rice <- do.call(rbind.data.frame,out)
post.rice.theta <- post.rice[,grepl('theta',colnames(post.rice))]
post.rice <- post.rice[,!grepl('theta',colnames(post.rice))]
agreement.rice <- agreementIndex(CRA=dat.rice$cra,CRAError=dat.rice$cra.error,calCurve='intcal20',theta=post.rice.theta)
# BPtoBCAD(rev(HPDinterval(mcmc(post.rice[,'a[1]']))))
# BPtoBCAD(rev(HPDinterval(mcmc(post.rice[,'a[2]']))))
save(rhats.rice,ess.rice,agreement.rice,post.rice,file=here('c14dates','04_ricearrival_results.RData'))
