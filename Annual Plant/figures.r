#calculate relationship between abundance and stabilization for all simulation results
#species 1 is the rare species
dat=read.csv("annplant_2spp.txt")
dat$ID=1:length(dat[,1])
dat$CoexistRank=NA
dat$CoexistRank[which(dat$medTC>0)]=0
dat$CoexistRank[which(dat$medTC>100)]=1
dat$CoexistRank[which(dat$medTC>1000)]=2
dat=dat[which(dat$l1>14),]
dat=dat[which(dat$Rare<=0.25),]

dat$S1cat=NA
dat$S1cat[which(dat$S1>=1)]="A"
dat$S1cat[which(dat$S1>=2)]="B"
dat$S1cat[which(dat$S1>=3)]="C"
dat$S1cat[which(dat$S1>=4)]="D"

dat$S2cat=NA
dat$S2cat[which(dat$S2>=1)]="A"
dat$S2cat[which(dat$S2>=2)]="B"
dat$S2cat[which(dat$S2>=3)]="C"
dat$S2cat[which(dat$S2>=4)]="D"

#effect of stabilization on persistence
plot(dat$S1,dat$medTC,log="y")
plot(dat$S2,dat$medTC,log="y")
t.test(S1~CoexistRank,data=dat[which(dat$CoexistRank!=2),])
t.test(S1~CoexistRank,data=dat[which(dat$CoexistRank!=1),])

t.test(dat$S1[which(dat$CoexistRank==0)],dat$S2[which(dat$CoexistRank==0)],paired=T)
t.test(dat$S1[which(dat$CoexistRank==1)],dat$S2[which(dat$CoexistRank==1)],paired=T)
t.test(dat$S1[which(dat$CoexistRank==2)],dat$S2[which(dat$CoexistRank==2)],paired=T)

#make correlation plots
plot(dat$S1[which(dat$cor<0)],dat$E1[which(dat$cor<0)],pch='-',cex=1.5, col=ifelse(dat$CoexistRank[which(dat$cor<0)]==0,"black",ifelse(dat$CoexistRank[which(dat$cor<0)]==1,"yellow","green")),ylim=c(0.5,1),xlab="Strength of stabilization",ylab="Fitness equivalence")
identify(dat$S1[which(dat$cor<0)],dat$E1[which(dat$cor<0)],labels=dat$ID[which(dat$cor<0)])
plot(dat$S1[which(dat$cor>=0)],dat$E1[which(dat$cor>=0)],pch='+',cex=1.5, col=ifelse(dat$CoexistRank[which(dat$cor>=0)]==0,"black",ifelse(dat$CoexistRank[which(dat$cor>=0)]==1,"yellow","green")),ylim=c(0.5,1),xlab="Strength of stabilization",ylab="Fitness equivalence")
identify(dat$S1[which(dat$cor>=0)],dat$E1[which(dat$cor>=0)],labels=dat$ID[which(dat$cor>=0)])

#make boxplot of coexistence vs strength of stabil by rank
boxplot(S1~CoexistRank,data=dat, names=c("1","1","1"), col=c("white","lightgrey","darkgrey"), xlab="Species", ylab="Strength of stabilization",at=c(1,4,7),xlim=c(0,9),ylim=c(1,13),range=0)
boxplot(S2~CoexistRank,data=dat, names=c("2","2","2"), col=c("white","lightgrey","darkgrey"), add=T,at=c(2,5,8),range=0)
abline(v=3,lwd=2)
abline(v=6,lwd=2)
text(1.2,12.5,labels="Persistence",cex=1)
text(1.4,12.1,labels="<100",cex=1)
text(4.5,12.5,labels="100<t<1000",cex=1)
text(7.4,12.5,labels=">1000",cex=1)

#make boxplot of coexistence vs diff in strength of stabil by rank
boxplot(Asy~CoexistRank,data=dat, names=c("<100","100<t<1000",">1000"),col=c("white","lightgrey","darkgrey"), xlab="Persistence time", ylab="Diff in strength of stabilization",range=0)

#plots of persistence vs stabilization
plot(dat$S1,dat$medTC,log="y",xlab="SOS", ylab="Median persistance time")

layout(matrix(c(1,2),1,2))
boxplot(medTC~S1cat,data=dat,xlab="Strength of stabilization",ylab="Median persistence time",log='y',boxwex=0.5)
boxplot(medTC~S2cat,data=dat,xlab="Strength of stabilization",ylab="Median persistence time",log='y',boxwex=0.5)

library(ggplot2)
source("C:\\Users\\gyenni\\Documents\\R\\multiplot.r")
means1=data.frame(S1cat=aggregate(dat$medTC,list(dat$S1cat),median)$Group.1,medTC=aggregate(dat$medTC,list(dat$S1cat),median)$x)
means2=data.frame(S2cat=aggregate(dat$medTC,list(dat$S2cat),median)$Group.1,medTC=aggregate(dat$medTC,list(dat$S2cat),median)$x)
#ggplot(dat, aes(medTC, ..density..)) + geom_histogram(binwidth=0.25,col="darkgrey",fill="darkgrey") + geom_density(aes(x = medTC),size=1) + geom_vline(aes(xintercept=medTC), data=means1, size=1) +
#facet_grid(S1cat ~ .) + scale_x_log() + theme_bw() + opts(strip.text.y = theme_text())

p1=ggplot(dat, aes(medTC, ..density..)) + geom_histogram(binwidth=0.5,col="darkgrey",fill="darkgrey") + geom_density(aes(x = medTC),size=1) +
  geom_vline(aes(xintercept=medTC), data=means1, size=1) + facet_wrap(~ S1cat,ncol=4) + scale_x_log(limits = c(3, 9200)) + theme_bw() + coord_flip() +
  xlab("median coexistence time") + opts(strip.text.x = theme_text(size=12,hjust=0),strip.background=theme_blank(), axis.title.y = theme_text(size = 14, angle = 90),axis.title.x = theme_text(size=14))

p2=ggplot(dat, aes(medTC, ..density..)) + geom_histogram(binwidth=0.5,col="darkgrey",fill="darkgrey") + geom_density(aes(x = medTC),size=1) +
  geom_vline(aes(xintercept=medTC), data=means2, size=1) + facet_wrap(~ S2cat,ncol=4) + scale_x_log(limits = c(3, 9200)) + theme_bw() + coord_flip() +
  xlab("median coexistence time") + opts( strip.text.x = theme_text(size = 12,hjust=0),strip.background=theme_blank(), axis.title.y = theme_text(size = 14, angle = 90),axis.title.x = theme_text(size=14))

multiplot(p1,p2,cols=1)

#simple plot for defense
library(plyr)
S1dat = ddply(dat,.(S1cat),summarise,med = median(medTC),medSE = sqrt(var(medTC)/length(medTC)))
S2dat = ddply(dat,.(S2cat),summarise,med = median(medTC),medSE = sqrt(var(medTC)/length(medTC)))
ggplot(dat) + ylab("median coexistence time") + xlab("Strength of NFD") +
  geom_pointrange(data=S1dat,aes(S1cat,med,ymin = med - medSE, ymax=med + medSE,colour=S1cat),size=2) +
  geom_pointrange(data=S2dat,aes(S2cat,med,ymin = med - medSE, ymax=med + medSE,colour=S2cat),size=2) +
  theme_bw() + geom_path(data=S1dat,aes(S1cat,med)) +
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
        axis.title.y = element_text(size = 18, angle = 90),axis.title.x = element_text(size=18))
  
#log-linear model of the effect of stabilization on persistence
summary(lm(log(medTC)~S1+E1+cor,data=dat))

coexist2=dat[which(dat$CoexistRank==2),]
coexist1=dat[which(dat$CoexistRank==1),]
extinct=dat[which(dat$CoexistRank==0),]
null=t.test(dat$Asy)
coex2=t.test(coexist2$Asy)
coex1=t.test(coexist1$Asy)
ext=t.test(extinct$Asy)
library(Hmisc)
errbar(c("<100","all","100<t<1000",">1000"),c(ext$estimate,null$estimate,coex1$estimate,coex2$estimate),c(ext$conf.int[2],null$conf.int[2],coex1$conf.int[2],coex2$conf.int[2]),c(ext$conf.int[1],null$conf.int[1],coex1$conf.int[1],coex2$conf.int[1]),ylab="Relationship between SOS and abundance")

#effect of l1, l2
plot(dat$l1,dat$l2,col=ifelse(dat$CoexistRank==0,"black",ifelse(dat$CoexistRank==1,"yellow","green")))
layout(matrix(c(1,2,3),1,3))
plot(dat$l1[which(dat$CoexistRank==2)]+runif(length(which(dat$CoexistRank==2))),dat$l2[which(dat$CoexistRank==2)]+runif(length(which(dat$CoexistRank==2))),xlim=c(10,21),ylim=c(10,21),col="black",pch=19,xlab="l1",ylab="l2",main="t>1000")
plot(dat$l1[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1))),dat$l2[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1))),col="grey",pch=19,xlab="l1",ylab="l2",main="100<t<1000")
plot(dat$l1[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0))),dat$l2[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0))),col="grey",xlab="l1",ylab="l2",main="t<100")

#deterministic data
detdat=read.csv("annplant_2spp_det_data.txt")
detrare=read.csv("annplant_2spp_det_rare.txt")

summary(glm(CoexistRank~S1+E1+cor,data=detrare,family=binomial(link="logit")))

dim(detdat)
dim(detrare)

dim(detrare)[1]/dim(detdat)[1]

length(which(detdat$Asy>=0))/length(detdat$Asy)
length(which(detrare$Asy>=0))/length(detrare$Asy)

length(which(detdat$Asy>=0))/length(detdat$Asy)
length(which(detrare$Asy>=0))/length(detrare$Asy)

t.test(detrare$S1[which(detrare$CoexistRank==1)],detrare$S2[which(detrare$CoexistRank==1)],paired=T)

t.test(S1~CoexistRank,data=detrare)
t.test(S2~CoexistRank,data=detrare)

t.test(Asy~CoexistRank,data=detdat)
t.test(Asy~CoexistRank,data=detrare)

#plot scenarios
  layout(matrix(c(1,2,3),1,3))
  plot(NA,NA,xlim=c(0,1), ylim=c(-1,1),xlab="",ylab="log growth rate",cex.lab=1.5,tcl=0.2, mgp=c(2,0.5,0))
  mtext(side=3,expression(paste("A)",nu > 0,sep=" ")),adj=0,line=0.5)
  abline(h=0,col="grey")
  abline(a=0.1,b=-1,lwd=2,lty=1)
  abline(a=0.6,b=-2,lwd=2,lty=2)
  abline(a=1.8,b=-3,lwd=2,lty=3)
  
  plot(NA,NA,xlim=c(0,1), ylim=c(-1,1),xlab="Frequency",ylab="",cex.lab=1.5,tcl=0.2, mgp=c(2,0.5,0))
  mtext(side=3,expression(paste("B)",nu %~~% 0,sep=" ")),adj=0,line=0.5)
  abline(h=0,col="grey")
  abline(a=0.2,b=-2,lwd=2,lty=1)
  abline(a=0.6,b=-2,lwd=2,lty=2)
  abline(a=1.2,b=-2,lwd=2,lty=3)
  
  plot(NA,NA,xlim=c(0,1), ylim=c(-1,1),xlab="",ylab="",cex.lab=1.5,tcl=0.2, mgp=c(2,0.5,0))
  mtext(side=3,expression(paste("C)",nu < 0,sep=" ")),adj=0,line=0.5)
  abline(h=0,col="grey")
  abline(a=0.3,b=-3,lwd=2,lty=1)
  abline(a=0.6,b=-2,lwd=2,lty=2)
  abline(a=0.6,b=-1,lwd=2,lty=3)

  #2-species cases
  layout(matrix(c(1,2,3),1,3))
  plot(NA,NA,xlim=c(0,1), ylim=c(-1,1),xlab="",ylab="log growth rate",cex.lab=1.5,tcl=0.2, mgp=c(2,0.5,0))
  mtext(side=3,expression(paste("A)",nu > 0,sep=" ")),adj=0,line=0.5)
  abline(h=0,col="grey")
  abline(a=0.2,b=-1,lwd=2,lty=1)
  abline(a=2.4,b=-3,lwd=2,lty=3)

  plot(NA,NA,xlim=c(0,1), ylim=c(-1,1),xlab="Frequency",ylab="",cex.lab=1.5,tcl=0.2, mgp=c(2,0.5,0))
  mtext(side=3,expression(paste("B)",nu %~~% 0,sep=" ")),adj=0,line=0.5)
  abline(h=0,col="grey")
  abline(a=0.4,b=-2,lwd=2,lty=1)
  abline(a=1.6,b=-2,lwd=2,lty=3)

  plot(NA,NA,xlim=c(0,1), ylim=c(-1,1),xlab="",ylab="",cex.lab=1.5,tcl=0.2, mgp=c(2,0.5,0))
  mtext(side=3,expression(paste("C)",nu < 0,sep=" ")),adj=0,line=0.5)
  abline(h=0,col="grey")
  abline(a=0.6,b=-3,lwd=2,lty=1)
  abline(a=0.8,b=-1,lwd=2,lty=3)

#plot simulation example

  layout(matrix(c(1,2),1))
  #S1>S2: #14385
  l1=dat$l1[7737]
  l2=dat$l2[7737]
  a11=dat$a11[7737]
  a12=dat$a12[7737]
  a21=dat$a21[7737]
  a22=dat$a22[7737]


  updateN = function(r_arg, Nself, N1, a_intra, a1){
    newN = (r_arg * Nself) /(1 + a_intra * Nself + a1 * N1 )
    newN=rpois(1,newN)         #demographic stochasticity
    return(newN)
  }
  # populations loop ---------------------------------
    N=data.frame("spp1"=5,"spp2"=5)  #initialize
    counter=1
    stopRun=F
    while(stopRun==F){
      counter=counter+1
    	N[counter,1] = updateN(l1,N[counter-1,1],N[counter-1,2],a11,a12)
    	N[counter,2] = updateN(l2,N[counter-1,2],N[counter-1,1],a22,a21)
      if(sum(N[counter,]==0)>=1) stopRun=T
      }
      dim(N)[1]
      matplot(N,type="l",lty=1:2,col=1,ylab="N",xlab="t")
      NL=N
      
  #S2>S1: #3982
  l1=dat$l1[1]
  l2=dat$l2[1]
  a11=dat$a11[1]
  a12=dat$a12[1]
  a21=dat$a21[1]
  a22=dat$a22[1]

  # populations loop ---------------------------------
    N=data.frame("spp1"=5,"spp2"=5)  #initialize
    counter=1
    stopRun=F
    while(stopRun==F){
      counter=counter+1
    	N[counter,1] = updateN(l1,N[counter-1,1],N[counter-1,2],a11,a12)
    	N[counter,2] = updateN(l2,N[counter-1,2],N[counter-1,1],a22,a21)
      if(sum(N[counter,]==0)>=1) stopRun=T
      }
      dim(N)[1]
      matplot(N,type="l",lty=1:2,col=1,ylab="N",xlab="t")
      NS=N
      
      
    matplot(NS[1:50,],type="l",lty=1:2,col=1,ylab="Density",xlab="",ylim=c(0,27),cex.lab=1.3,tcl=0.2, mgp=c(2,0.5,0))
    mtext(side=3,"A",adj=0,line=0.5)
    matplot(NL[1:50,],type="l",lty=1:2,col=1,ylab="",xlab="",ylim=c(0,27),cex.lab=1.3,tcl=0.2, mgp=c(2,0.5,0))
    mtext(side=3,"B",adj=0,line=0.5)
    par(xpd=NA)
    legend(locator(1),legend="Year",bty='n',cex=1.3)
    