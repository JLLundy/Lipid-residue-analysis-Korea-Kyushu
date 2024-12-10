# Load libraries and Read Data ----
library(rcarbon)
library(dplyr)
library(nimbleCarbon)
library(coda)
library(parallel)
library(here)

# Arrval dates
load(here('arrivalestimates','ricearrival_res.RData'))
load(here('arrivalestimates','milletsarrival_res.RData'))

# Read and Calibrate Dates
korea.dates.rice <- read.csv(here('demography','S0033822217001229sup001.csv')) |> rename(c14age=date,c14error=error)
korea.caldates.rice <- calibrate(korea.dates.rice$c14age,korea.dates.rice$c14error)
load(url('https://github.com/ercrema/NeolithicKoreaDemography/raw/master/R_image_files/koreanC14.RData'))
korea.dates.millets <- koreaC14
korea.caldates.millets <- calibrate(korea.dates.millets$c14age,korea.dates.millets$c14error)


# Extract median estimated arrival time for rice
(ricearrival.sk <- median(post.rice[,'a[2]']) |> round())
# 3304
(milletsarrival.sk <- median(post.millets[,'a[2]']) |> round())
# 6116

# Subset dates with probability mass over 0.5 after the arrival dates
i500.rice <- which.CalDates(korea.caldates.rice,BP < ricearrival.sk & BP > (ricearrival.sk-500),p=0.5)
i500.millets <- which.CalDates(korea.caldates.millets,BP < milletsarrival.sk & BP > (milletsarrival.sk-500),p=0.5)

korea.dates.500.rice <- korea.dates.rice[i500.rice,]
korea.caldates.500.rice <- korea.caldates.rice[i500.rice]
korea.dates.500.rice$med <- medCal(korea.caldates.500.rice)
korea.dates.500.rice$SiteID <- as.numeric(as.factor(korea.dates.500.rice$Site))
korea.bins.500.rice <- binPrep(sites=korea.dates.500.rice$SiteID,korea.caldates.500.rice,h=50)
thin500.rice <- thinDates(ages=korea.dates.500.rice$c14age,errors=korea.dates.500.rice$c14error,bins=korea.bins.500.rice,size=1,thresh=1,seed=123,method='splitsample')
korea.dates.500.thinned.rice <- korea.dates.500.rice[thin500.rice,]

korea.dates.500.millets <- korea.dates.millets[i500.millets,]
korea.caldates.500.millets <- korea.caldates.millets[i500.millets]
korea.dates.500.millets$med <- medCal(korea.caldates.500.millets)
korea.dates.500.millets$SiteID <- as.numeric(as.factor(korea.dates.500.millets$site_id))
korea.bins.500.millets <- binPrep(sites=korea.dates.500.millets$SiteID,korea.caldates.500.millets,h=50)
thin500.millets <- thinDates(ages=korea.dates.500.millets$c14age,errors=korea.dates.500.millets$c14error,bins=korea.bins.500.millets,size=1,thresh=1,seed=123,method='splitsample')
korea.dates.500.thinned.millets <- korea.dates.500.millets[thin500.millets,]


# Define constants and data ----
# data
dat.rice  <- dat.millets  <- list()
dat.rice$cra <- korea.dates.500.thinned.rice$c14age
dat.rice$cra.error <- korea.dates.500.thinned.rice$c14error
dat.millets$cra <- korea.dates.500.thinned.millets$c14age
dat.millets$cra.error <- korea.dates.500.thinned.millets$c14error


# constants
data(intcal20)
constants.base <- list()
constants.base$calBP  <- intcal20$CalBP
constants.base$C14BP <- intcal20$C14Age
constants.base$C14Err <- intcal20$C14Age.sigma
constants.rice <- constants.millets <- constants.base
constants.rice$n <- nrow(korea.dates.500.thinned.rice)
constants.rice$a <- ricearrival.sk
constants.rice$b <- ricearrival.sk - 500
constants.millets$n <- nrow(korea.dates.500.thinned.millets)
constants.millets$a <- milletsarrival.sk
constants.millets$b <- milletsarrival.sk - 500

#theta init
theta.rice <- korea.dates.500.thinned.rice$med
if(any(theta.rice > constants.rice$a | theta.rice < constants.rice$b))
{
	theta.rice[which(theta,rice>constants.rice$a)] <- constants.rice$a - 1
	theta.rice[which(theta,rice<constants.rice$b)] <- constants.rice$b + 1
}

theta.millets <- korea.dates.500.thinned.millets$med
if(any(theta.millets > constants.millets$a | theta.millets < constants.millets$b))
{
	theta.millets[which(theta,millets>constants.millets$a)] <- constants.millets$a - 1
	theta.millets[which(theta,millets<constants.millets$b)] <- constants.millets$b + 1
}

# Define run function ----
runFun <- function(seed,d,constants,theta.inits,niter,nburnin,thin)
{
	library(nimbleCarbon)
	exponential <- nimbleCode({
		for (i in 1:n){
			theta[i] ~ dExponentialGrowth(a=a,b=b,r=r);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14Err[]);
			sd[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}
		r ~ dnorm(mean=0,sd=0.001)
	})
	set.seed(seed)
	inits <- list()
	inits$theta=theta.inits
	inits$r  <- rnorm(1,mean=0,sd=0.001)
	#MCMC
	model <- nimbleModel(exponential, constants = constants, data = d, inits = inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model, monitors = c('r','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC, project = cModel)
	samples <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(samples)
}

# Fit Models ----
niter  <- 500000
nburnin <- 250000
thin <- 25
seeds  <- 1:4
ncores <- 4

cl <- makeCluster(ncores)
out.rice  <-  parLapply(cl=cl,X=seeds,fun=runFun,d=dat.rice,constants=constants.rice,theta.inits=theta.rice,niter=niter,nburnin=nburnin,thin=thin)
rhats.rice <- gelman.diag(out.rice)
any(rhats.rice[[1]][,1]>1.01)
post.rice <- do.call(rbind.data.frame,out.rice)
post.rice.r <- post.rice$r
stopCluster(cl)

cl <- makeCluster(ncores)
out.millets  <-  parLapply(cl=cl,X=seeds,fun=runFun,d=dat.millets,constants=constants.millets,theta.inits=theta.millets,niter=niter,nburnin=nburnin,thin=thin)
rhats.millets <- gelman.diag(out.millets)
any(rhats.millets[[1]][,1]>1.01)
post.millets <- do.call(rbind.data.frame,out.millets)
post.millets.r <- post.millets$r
stopCluster(cl)

save(post.rice.r,rhats.rice,post.millets.r,rhats.millets,file=here('results','growth_res_korea.RDataa'))

