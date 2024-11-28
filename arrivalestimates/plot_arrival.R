# Load Library and Posteriors ----
library(here)
library(coda)

load(here('scripts','ricearrival_res.RData'))
load(here('scripts','milletsarrival_res.RData'))
load(here('scripts','setariaarrival_res.RData'))
load(here('scripts','panicumarrival_res.RData'))

# Utility function ----
posterior.bar <- function(x,y,width=0.4,prob1=0.5,prob2=0.9,fill='grey',alpha=0.5)
{
	x.lo1 <- HPDinterval(mcmc(x),prob=prob1)[1]
	x.hi1 <- HPDinterval(mcmc(x),prob=prob1)[2]
	x.lo2 <- HPDinterval(mcmc(x),prob=prob2)[1]
	x.hi2 <- HPDinterval(mcmc(x),prob=prob2)[2]
	med <- median(x)
	rect(xleft=x.lo1,xright=x.hi1,ybottom=y-width/2,ytop=y+width/2,border=NA,col=fill)
	rect(xleft=x.lo2,xright=x.lo1,ybottom=y-width/2,ytop=y+width/2,border=NA,col=adjustcolor(fill,alpha))
	rect(xleft=x.hi1,xright=x.hi2,ybottom=y-width/2,ytop=y+width/2,border=NA,col=adjustcolor(fill,alpha))
# 	points(x=med,y=y,pch=20,cex=cex)
}


posterior.plot <- function(x,prob=0.9,bw='SJ',hpd.col='lightblue',line.col='darkgrey',add=FALSE,...)
{
	x.lo <- HPDinterval(mcmc(x),prob=prob)[1]
	x.hi <- HPDinterval(mcmc(x),prob=prob)[2]
	dens <- density(x,bw=bw)
	if (add==FALSE){plot(dens$x,dens$y,type='n',...)}
	hpdi.x <- dens$x[which(dens$x>=x.lo&dens$x<=x.hi)]
	hpdi.y <- dens$y[which(dens$x>=x.lo&dens$x<=x.hi)]
	polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col=hpd.col)
	lines(dens,col=line.col)
}

posterior.bar <- function(x,prob=0.9,calendar='BCAD',yp=0.8,ys=0.05)
{
	require(rcarbon)
	ygap  <- c(par('usr')[4]-par('usr')[3])*yp
	ygap2  <- c(par('usr')[4]-par('usr')[3])*(yp+ys)
	x.lo <- HPDinterval(mcmc(x),prob=prob)[1]
	x.hi <- HPDinterval(mcmc(x),prob=prob)[2]
	lines(c(x.lo,x.hi),c(ygap,ygap))
	points(c(x.lo,x.hi),c(ygap,ygap),pch=20)
	labs <- c(x.hi,x.lo) |> round()
	if(calendar=='BCAD'){labs <- abs(BPtoBCAD(labs))}
	text(c(x.hi,x.lo),c(ygap2,ygap2),labels=labs)
}


col1='#2081f9'
col2='#f99820'
width=0.3
pdf('test.pdf',width=8,height=3)
layout(matrix(1:2,nrow=1,ncol=2),width=c(0.6,0.4))
par(mar=c(3,0,0,0))
plot(NULL,xlim=c(8000,2500),ylim=c(0.5,4.5),xlab='',ylab='',axes=F)
posterior.bar(x=post.rice[,'a[1]'],y=4.2,fill=col1,width=width)
posterior.bar(x=post.rice[,'a[2]'],y=3.8,fill=col2,width=width)
posterior.bar(x=post.millets[,'a[1]'],y=3.2,fill=col1,width=width)
posterior.bar(x=post.millets[,'a[2]'],y=2.8,fill=col2,width=width)
posterior.bar(x=post.setaria[,'a[1]'],y=2.2,fill=col1,width=width)
posterior.bar(x=post.setaria[,'a[2]'],y=1.8,fill=col2,width=width)
posterior.bar(x=post.panicum[,'a[1]'],y=1.2,fill=col1,width=width)
posterior.bar(x=post.panicum[,'a[2]'],y=0.8,fill=col2,width=width)
abline(h=c(1.5,2.5,3.5),col='lightgrey')
text(x=BCADtoBP(-5200),y=4.2,labels='Rice',cex=0.8)
text(x=BCADtoBP(-5200),y=3.2,labels='Millets (combined)',cex=0.8)
text(x=BCADtoBP(-5200),y=2.2,labels='Foxtail Millet',cex=0.8)
text(x=BCADtoBP(-5200),y=1.2,labels='Broomcorn Millet',cex=0.8)
axis(1,at=BCADtoBP(seq(-6000,-500,500)),labels=abs(seq(-6000,-500,500)),cex.axis=0.6,mgp=c(1.6,0.1,0),tck=-0.03)
axis(1,at=BCADtoBP(seq(-6000,-500,100)),labels=NA,tck=-0.01)
mtext(side=1,text='BCE',line=1,cex=0.6)
par(mar=c(3,0,0,0))
plot(NULL,xlim=c(0,6500),ylim=c(0.5,4.5),xlab='',ylab='',axes=F)
posterior.bar(x=post.rice[,'a[2]']-post.rice[,'a[1]'],y=4,fill='grey',width=width)
posterior.bar(x=post.millets[,'a[2]']-post.millets[,'a[1]'],y=3,fill='grey',width=width)
posterior.bar(x=post.setaria[,'a[2]']-post.setaria[,'a[1]'],y=2,fill='grey',width=width)
posterior.bar(x=post.panicum[,'a[2]']-post.panicum[,'a[1]'],y=1,fill='grey',width=width)
axis(1,cex.axis=0.6,mgp=c(1.6,0.1,0),tck=-0.03)
axis(1,at=seq(0,6500,200),tck=-0.01,labels=NA)
abline(h=c(1.5,2.5,3.5),col='lightgrey')
mtext(side=1,tex='Delay (years)',line=1,cex=0.6)
legend('topright',legend=c('South Korea','Japan'),fill=c(col2,col1),cex=0.8,bty='n')
dev.off()


