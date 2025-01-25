library(here)
library(parallel)
library(coda)
library(rcarbon)
library(nimbleCarbon)
load(here('c14dates','02_millets_run.RData'))

#MCMC settings ----
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12,34,56,78)
niter  <- 2000000
nburnin  <- 1000000
thin  <- 100

# Run Model
out  <-  parLapply(cl = cl, X = seeds, fun = runscript, d = dat.millets, constants = constants.millets, inits = inits.millets, niter = niter, nburnin = nburnin,thin = thin)
stop(cl)

rhats.millets <- gelman.diag(out,multivariate=T)
rhats.millets[[1]][which(rhats.millets[[1]][,1]>1.01),]
ess.millets <- effectiveSize(out)
post.millets <- do.call(rbind.data.frame,out)
post.millets.theta <- post.millets[,grepl('theta',colnames(post.millets))]
post.millets <- post.millets[,!grepl('theta',colnames(post.millets))]
agreement.millets <- agreementIndex(CRA=dat.millets$cra,CRAError=dat.millets$cra.error,calCurve='intcal20',theta=post.millets.theta)
save(rhats.millets,ess.millets,agreement.millets,post.millets,file=here('c14dates','04_milletsarrival_results.RData'))



