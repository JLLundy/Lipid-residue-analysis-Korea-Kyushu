library(here)
library(parallel)
library(coda)
library(rcarbon)
load(here('arrivalestimates','setaria_run.RData'))

#MCMC settings ----
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12,34,56,78)
niter  <- 6000000
nburnin  <- 3000000
thin  <- 300

# Run Model
out  <-  parLapply(cl = cl, X = seeds, fun = runscript, d = dat.setaria, constants = constants.setaria, inits = inits.setaria, niter = niter, nburnin = nburnin,thin = thin)
stop(cl)

rhats.setaria <- gelman.diag(out,multivariate=T)
ess.setaria <- effectiveSize(out)
post.setaria <- do.call(rbind.data.frame,out)
BPtoBCAD(rev(HPDinterval(mcmc(post.setaria[,'a[1]']))))
BPtoBCAD(rev(HPDinterval(mcmc(post.setaria[,'a[2]']))))
agreement.setaria  <- agreementIndex(CRA=dat.setaria$cra,CRAError=dat.setaria$cra.error,calCurve='intcal20',theta=post.setaria.theta)
save(rhats.setaria,ess.setaria,agreement.setaria,post.setaria,file=here('arrivalestimates','setariaarrival_res.RData'))


