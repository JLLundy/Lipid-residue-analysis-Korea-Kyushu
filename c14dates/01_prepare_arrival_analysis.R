# Load Libraries and Read Files ----
library(here)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
d  <-  read.csv(here('data','S5.csv')) |> subset(!exclude=='YES')
cal.dates <- calibrate(d$CRA,d$CRAError)
d$med <- medCal(cal.dates)
d <- d[which.CalDates(cal.dates,BP>1500,p=0.5),]
data(intcal20)

d.rice <- subset(d,Material=='rice') |> select(SiteName,Region,CRA,CRAError,med)
d.rice$SiteID  <- as.numeric(as.factor(d.rice$SiteName))
d.rice$RegionID  <- as.numeric(as.factor(d.rice$Region))

d.millets <- subset(d,Material%in%c('foxtail millet','broomcorn millet')) |> select(SiteName,Region,CRA,CRAError,med)
d.millets$SiteID  <- as.numeric(as.factor(d.millets$SiteName))
d.millets$RegionID  <- as.numeric(as.factor(d.millets$Region))

# Setup Data and Constants ----
# data
dat.rice <- dat.millets <- list()

dat.rice$cra <- d.rice$CRA
dat.rice$cra.error <- d.rice$CRAError

dat.millets$cra <- d.millets$CRA
dat.millets$cra.error <- d.millets$CRAError

dat.rice$constraint.region <- dat.millets$constraint.region <- c(1,1)

# Constants and inits
constants.rice  <- constants.millets <- list()
inits.rice  <- inits.millets <- list()

constants.rice$n <- nrow(d.rice)
constants.rice$n.region <- 2
constants.rice$n.sites <- length(unique(d.rice$SiteID))
constants.rice$site.id <- d.rice$SiteID
constants.rice$region.id <- select(d.rice,SiteID,RegionID) |> unique() |> arrange(SiteID) |> select(RegionID)|>as.matrix()|>as.numeric()
inits.rice$theta <- d.rice$med
inits.rice$alpha  <- select(d.rice,SiteID,med) |> unique() |> aggregate(med~SiteID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 50
inits.rice$delta  <- select(d.rice,SiteID,med) |> unique() |> aggregate(med~SiteID,function(x){abs(min(x)-max(x))},data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.rice$a <- select(d.rice,RegionID,med) |> unique() |> aggregate(med~RegionID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.rice$b <- select(d.rice,RegionID,med) |> unique() |> aggregate(med~RegionID,min,data=_) |> select(med) |> as.matrix() |> as.numeric() - 100

constants.millets$n <- nrow(d.millets)
constants.millets$n.region <- 2
constants.millets$n.sites <- length(unique(d.millets$SiteID))
constants.millets$site.id <- d.millets$SiteID
constants.millets$region.id <- select(d.millets,SiteID,RegionID) |> unique() |> arrange(SiteID) |> select(RegionID)|>as.matrix()|>as.numeric()
inits.millets$theta <- d.millets$med
inits.millets$alpha  <- select(d.millets,SiteID,med) |> unique() |> aggregate(med~SiteID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 50
inits.millets$delta  <- select(d.millets,SiteID,med) |> unique() |> aggregate(med~SiteID,function(x){abs(min(x)-max(x))},data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.millets$a <- select(d.millets,RegionID,med) |> unique() |> aggregate(med~RegionID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.millets$b <- select(d.millets,RegionID,med) |> unique() |> aggregate(med~RegionID,min,data=_) |> select(med) |> as.matrix() |> as.numeric() - 100

constants.rice$calBP <- constants.millets$calBP <- intcal20$CalBP
constants.rice$C14BP <- constants.millets$C14BP <- intcal20$C14Age
constants.rice$C14err <- constants.millets$C14err  <- intcal20$C14Age.sigma

# Define Runscript ----

runscript <- function(seed,d,inits,constants,nburnin,thin,niter)
{
	library(nimbleCarbon)
	model <- nimbleCode({
		for (k in 1:n.sites)
		{
			delta[k] ~ dgamma(gamma1,(gamma1-1)/gamma2)
			alpha[k] ~ dunif(max=a[region.id[k]],min=b[region.id[k]])
		}

		for (i in 1:n){
			theta[i] ~ dunif(alpha[site.id[i]]-delta[site.id[i]]+1,alpha[site.id[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sd[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=sd[i]);
		}

		#Prior for each region
		for (j in 1:n.region){
			a[j] ~ dunif(50,10000)
			b[j] ~ dunif(50,10000)
			constraint.region[j] ~ dconstraint(a[j] > b[j])
		}

		# Hyperprior for duration
		gamma1 ~ dunif(1,20)
		gamma2 ~ T(dnorm(mean=200,sd=100),1,500)
	})

	inits$gamma1 <- 10
	inits$gamma2 <- 200

	# Compile and Run model	
	model <- nimbleModel(model,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,control=list(adaptInterval=20000,adaptFactorExponent=0.1))
	conf$addMonitors(c('theta','delta','alpha'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
}

save(dat.rice,constants.rice,inits.rice,runscript,file=here('c14dates','02_rice_data.RData'))
save(dat.millets,constants.millets,inits.millets,runscript,file=here('c14dates','02_millets_run.RData'))
