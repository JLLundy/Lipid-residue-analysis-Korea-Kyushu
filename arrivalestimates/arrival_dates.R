# Load Libraries and Read Files ----
library(here)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
d  <-  read.csv(here('arrivalestimates','combined_final.csv')) |> subset(!Exclude=='YES')
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

d.setaria <- subset(d,Material=='foxtail millet') |> select(SiteName,Region,CRA,CRAError,med)
d.setaria$SiteID  <- as.numeric(as.factor(d.setaria$SiteName))
d.setaria$RegionID  <- as.numeric(as.factor(d.setaria$Region))

d.panicum <- subset(d,Material=='broomcorn millet') |> select(SiteName,Region,CRA,CRAError,med)
d.panicum$SiteID  <- as.numeric(as.factor(d.panicum$SiteName))
d.panicum$RegionID  <- as.numeric(as.factor(d.panicum$Region))

# Setup Data and Constants ----
# data
dat.rice <- dat.millets <- dat.setaria <- dat.panicum <- list()

dat.rice$cra <- d.rice$CRA
dat.rice$cra.error <- d.rice$CRAError

dat.millets$cra <- d.millets$CRA
dat.millets$cra.error <- d.millets$CRAError

dat.setaria$cra <- d.setaria$CRA
dat.setaria$cra.error <- d.setaria$CRAError

dat.panicum$cra <- d.panicum$CRA
dat.panicum$cra.error <- d.panicum$CRAError

dat.rice$constraint.region <- dat.millets$constraint.region <- dat.setaria$constraint.region <- dat.panicum$constraint.region <- c(1,1)

# Constants and inits
constants.rice  <- constants.millets <- constants.setaria <- constants.panicum <- list()
inits.rice  <- inits.millets <- inits.setaria <- inits.panicum <- list()

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

constants.setaria$n <- nrow(d.setaria)
constants.setaria$n.region <- 2
constants.setaria$n.sites <- length(unique(d.setaria$SiteID))
constants.setaria$site.id <- d.setaria$SiteID
constants.setaria$region.id <- select(d.setaria,SiteID,RegionID) |> unique() |> arrange(SiteID) |> select(RegionID)|>as.matrix()|>as.numeric()
inits.setaria$theta <- d.setaria$med
inits.setaria$alpha  <- select(d.setaria,SiteID,med) |> unique() |> aggregate(med~SiteID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 50
inits.setaria$delta  <- select(d.setaria,SiteID,med) |> unique() |> aggregate(med~SiteID,function(x){abs(min(x)-max(x))},data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.setaria$a <- select(d.setaria,RegionID,med) |> unique() |> aggregate(med~RegionID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.setaria$b <- select(d.setaria,RegionID,med) |> unique() |> aggregate(med~RegionID,min,data=_) |> select(med) |> as.matrix() |> as.numeric() - 100

constants.panicum$n <- nrow(d.panicum)
constants.panicum$n.region <- 2
constants.panicum$n.sites <- length(unique(d.panicum$SiteID))
constants.panicum$site.id <- d.panicum$SiteID
constants.panicum$region.id <- select(d.panicum,SiteID,RegionID) |> unique() |> arrange(SiteID) |> select(RegionID)|>as.matrix()|>as.numeric()
inits.panicum$theta <- d.panicum$med
inits.panicum$alpha  <- select(d.panicum,SiteID,med) |> unique() |> aggregate(med~SiteID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 50
inits.panicum$delta  <- select(d.panicum,SiteID,med) |> unique() |> aggregate(med~SiteID,function(x){abs(min(x)-max(x))},data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.panicum$a <- select(d.panicum,RegionID,med) |> unique() |> aggregate(med~RegionID,max,data=_) |> select(med) |> as.matrix() |> as.numeric() + 100
inits.panicum$b <- select(d.panicum,RegionID,med) |> unique() |> aggregate(med~RegionID,min,data=_) |> select(med) |> as.matrix() |> as.numeric() - 100

constants.rice$calBP <- constants.millets$calBP <- constants.setaria$calBP <- constants.panicum$calBP <- intcal20$CalBP
constants.rice$C14BP <- constants.millets$C14BP <- constants.setaria$C14BP <- constants.panicum$C14BP <- intcal20$C14Age
constants.rice$C14err <- constants.millets$C14err <- constants.setaria$C14err <- constants.panicum$C14err <- intcal20$C14Age.sigma

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
	conf <- configureMCMC(model)
	conf$addMonitors(c('theta','delta','alpha'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed=seed) 
}

save(dat.rice,constants.rice,inits.rice,runscript,file=here('arrivaldates','rice_run.RData'))
save(dat.millets,constants.millets,inits.millets,runscript,file=here('arrivaldates','millets_run.RData'))
save(dat.setaria,constants.setaria,inits.setaria,runscript,file=here('arrivaldates','setaria_run.RData'))
save(dat.panicum,constants.panicum,inits.panicum,runscript,file=here('arrivaldates','panicum_run.RData'))


