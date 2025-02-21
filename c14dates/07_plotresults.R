# Load Library and Posteriors ----
library(here)
library(coda)
library(rcarbon)

# Demography estimates ----

# Estimates from Crema et al 2024 (N.Kyushu) doi:10.15184/aqy.2024.146
load(url('https://github.com/ercrema/yayoi_demo/raw/refs/heads/main/results/icar_c14doubleRes500.RData'))
icar.res.c14.500  <-  icar.c14double500[,1:16] * 100
res500.90  <- apply(icar.res.c14.500,2,function(x){HPDinterval(as.mcmc(x),0.90)})
res500.50  <- apply(icar.res.c14.500,2,function(x){HPDinterval(as.mcmc(x),0.50)})
res500.med  <- apply(icar.res.c14.500,2,median)
nkyushu.growth.90 <- res500.90[,'s2[1]']
nkyushu.growth.50 <- res500.50[,'s2[1]']
nkyushu.growth.med <- res500.med['s2[1]']
start.nkyushu <- BCADtoBP(-1039)
end.nkyushu <- start.nkyushu - 500

# Estimates for Korea
load(here('c14dates','06_koreagrowth_results.RData'))
post.millets.r.hpdi.90  <- HPDinterval(as.mcmc(post.millets.r),0.9)
post.millets.r.hpdi.50  <- HPDinterval(as.mcmc(post.millets.r),0.5)
post.millets.r.hpdi.med  <- median(post.millets.r)
post.rice.r.hpdi.90  <- HPDinterval(as.mcmc(post.rice.r),0.9)
post.rice.r.hpdi.50  <- HPDinterval(as.mcmc(post.rice.r),0.5)
post.rice.r.hpdi.med  <- median(post.rice.r)


# Arrival Estimates ----
d <- read.csv(here('data','S5.csv')) |> subset(!exclude=='YES')
d$med <- medCal(calibrate(d$CRA,d$CRAError))
d.korea.rice  <- subset(d,Region=='South Korea'&Material=='rice')
d.japan.rice <- subset(d,Region=='Japan'&Material=='rice')
d.korea.millets  <- subset(d,Region=='South Korea'&Material!='rice')
d.japan.millets <- subset(d,Region=='Japan'&Material!='rice')


load(here('c14dates','04_milletsarrival_results.RData'))
load(here('c14dates','04_ricearrival_results.RData'))

# Utility function ----
posterior.bar <- function(x,y,width=0.4,prob1=0.5,prob2=0.9,fill='grey',lcol='black',alpha=0.5)
{
	x.lo1 <- HPDinterval(mcmc(x),prob=prob1)[1]
	x.hi1 <- HPDinterval(mcmc(x),prob=prob1)[2]
	x.lo2 <- HPDinterval(mcmc(x),prob=prob2)[1]
	x.hi2 <- HPDinterval(mcmc(x),prob=prob2)[2]
	med <- median(x)
	rect(xleft=x.lo1,xright=x.hi1,ybottom=y-width/2,ytop=y+width/2,border=NA,col=fill)
	rect(xleft=x.lo2,xright=x.lo1,ybottom=y-width/2,ytop=y+width/2,border=NA,col=adjustcolor(fill,alpha))
	rect(xleft=x.hi1,xright=x.hi2,ybottom=y-width/2,ytop=y+width/2,border=NA,col=adjustcolor(fill,alpha))
	lines(x=c(med,med),y=c(y-width/2,y+width/2),lwd=2,col=lcol)
}


# Estimates for in-text mentions

# Arrival time rice in Japan
HPDinterval(mcmc(post.rice[,'a[1]'])) |> BPtoBCAD()
# Arrival time rice in Korea
HPDinterval(mcmc(post.rice[,'a[2]'])) |> BPtoBCAD()
# Difference arrival time in rice
HPDinterval(mcmc(post.rice[,'a[2]']-post.rice[,'a[1]']))

# Arrival time millets in Japan
HPDinterval(mcmc(post.millets[,'a[1]'])) |> BPtoBCAD()
# Arrival time millets in Korea
HPDinterval(mcmc(post.millets[,'a[2]'])) |> BPtoBCAD()
# Difference arrival time in millets
HPDinterval(mcmc(post.millets[,'a[2]']-post.millets[,'a[1]']))


col1='#2081f9'
col2='#f99820'
width=0.03
pdf(here('figures','c14_analysis.pdf'),width=8,height=5)
# layout(matrix(1:2,nrow=1,ncol=2),width=c(0.6,0.4))
par(mar=c(3,4,0,2),xpd=T)
plot(NULL,xlim=c(7500,2500),ylim=c(-0.27,0.55),xlab='',ylab='',axes=F)

rect(xleft=median(post.rice[,'a[2]']),xright=median(post.rice[,'a[2]'])-500,ybottom=post.rice.r.hpdi.90[1]*100,ytop=post.rice.r.hpdi.90[2]*100,border=NA,col=adjustcolor(col2,0.5))
rect(xleft=median(post.rice[,'a[2]']),xright=median(post.rice[,'a[2]'])-500,ybottom=post.rice.r.hpdi.50[1]*100,ytop=post.rice.r.hpdi.50[2]*100,border=NA,col=adjustcolor(col2,0.9))
lines(x=c(median(post.rice[,'a[2]']),median(post.rice[,'a[2]'])-500),y=c(post.rice.r.hpdi.med,post.rice.r.hpdi.med)*100,lwd=4)
rect(xleft=BCADtoBP(-1039),xright=BCADtoBP(-1039+500),ybottom=nkyushu.growth.90[1],ytop=nkyushu.growth.90[2],border=NA,col=adjustcolor(col1,0.5))
rect(xleft=BCADtoBP(-1039),xright=BCADtoBP(-1039+500),ybottom=nkyushu.growth.50[1],ytop=nkyushu.growth.50[2],border=NA,col=adjustcolor(col1,0.9))
lines(x=c(BCADtoBP(c(-1039,-1039+500))),y=c(nkyushu.growth.med,nkyushu.growth.med),lwd=4)

rect(xleft=median(post.millets[,'a[2]']),xright=median(post.millets[,'a[2]'])-500,ybottom=post.millets.r.hpdi.90[1]*100,ytop=post.millets.r.hpdi.90[2]*100,border=NA,col=adjustcolor(col2,0.5))
rect(xleft=median(post.millets[,'a[2]']),xright=median(post.millets[,'a[2]'])-500,ybottom=post.millets.r.hpdi.50[1]*100,ytop=post.millets.r.hpdi.50[2]*100,border=NA,col=adjustcolor(col2,0.9))
lines(x=c(median(post.millets[,'a[2]']),median(post.millets[,'a[2]'])-500),y=c(post.millets.r.hpdi.med,post.millets.r.hpdi.med)*100,lwd=4)
lines(x=c(median(post.rice[,'a[2]']),median(post.rice[,'a[2]'])-500),y=c(post.rice.r.hpdi.med,post.rice.r.hpdi.med)*100,lwd=4)


posterior.bar(x=post.rice[,'a[1]'],y=-0.1,fill=col1,width=width)
points(d.japan.rice$med,y=rep(-0.1,nrow(d.japan.rice)),pch=18,col=adjustcolor('black',0.2))
posterior.bar(x=post.rice[,'a[2]'],y=-0.15,fill=col2,width=width)
points(d.korea.rice$med,y=rep(-0.15,nrow(d.korea.rice)),pch=18,col=adjustcolor('black',0.2))
posterior.bar(x=post.millets[,'a[1]'],y=-0.2,fill=col1,width=width)
points(d.japan.millets$med,y=rep(-0.2,nrow(d.japan.millets)),pch=18,col=adjustcolor('black',0.2))
posterior.bar(x=post.millets[,'a[2]'],y=-0.25,fill=col2,width=width)
points(d.korea.millets$med,y=rep(-0.25,nrow(d.korea.millets)),pch=18,col=adjustcolor('black',0.2))

abline(h=c(-0.075),col='darkgrey')
abline(h=c(-0.175),col='lightgrey',lty=2)
text(x=8100,y=-0.13,labels='Rice',cex=0.8,srt=90)
text(x=8100,y=-0.235,labels='Millets',cex=0.8,srt=90)
axis(1,at=BCADtoBP(seq(-6000,-500,500)),labels=abs(seq(-6000,-500,500)),cex.axis=0.8,mgp=c(1.6,0.1,0),tck=-0.01)
axis(1,at=BCADtoBP(seq(-6000,-500,100)),labels=NA,tck=-0.007)
mtext(side=1,text='BCE',line=1,cex=0.8)
axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,0.1,0.2,0.3,0.4,0.5),cex.axis=0.8,mgp=c(1.5,0.4,0),tck=-0.01,las=2)
text(x=8100,y=0.25,'Annual growth rate (%)',cex=0.9,srt=90)
legend('topleft',legend=c('Korean Peninsula','Japan'),fill=c(col2,col1),cex=1,bty='n')
dev.off()

