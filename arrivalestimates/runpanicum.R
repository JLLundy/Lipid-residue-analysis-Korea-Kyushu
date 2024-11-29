library(here)
library(parallel)
library(coda)
library(rcabon)

load(here('panicum_run.RData'))

#MCMC settings ----
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12,34,56,78)
niter  <- 6000000
nburnin  <- 3000000
thin  <- 300

# Run Model
out  <-  parLapply(cl = cl, X = seeds, fun = runscript, d = dat.panicum, constants = constants.panicum, inits = inits.panicum, niter = niter, nburnin = nburnin,thin = thin)
stop(cl)

rhats.panicum <- gelman.diag(out,multivariate=T)
ess.panicum <- effectiveSize(out)
post.panicum <- do.call(rbind.data.frame,out)
BPtoBCAD(rev(HPDinterval(mcmc(post.panicum[,'a[1]']))))
BPtoBCAD(rev(HPDinterval(mcmc(post.panicum[,'a[2]']))))
agreement.panicum  <- agreementIndex(CRA=dat.panicum$cra,CRAError=dat.panicum$cra.error,calCurve='intcal20',theta=post.panicum.theta)
save(rhats.panicum,ess.panicum,agreement.panicum,post.panicum,file=here('arrivalestimates','panicumarrival_res.RData'))





