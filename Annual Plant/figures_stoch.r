#calculate relationship between abundance and stabilization for all simulation results
#species 1 is the inferior species
dat=read.csv("annplant_2spp_stoch_rare.txt")

#dat=dat[which(dat$Rank1==2),]
#dat=dat[which(dat$l1>10),]

#calculate correlation btw abund and stabilization
corrs=rep(NA,length(dat[,1]))
for(i in 1:length(dat[,1])) {
corrs[i]=cov(x=c(dat[i,6],dat[i,8]),y=c(dat[i,21],dat[i,23]))  }
dat$cor=corrs

#effect of stabilization on persistence
plot(dat$S1,dat$medTC,log="y")
plot(dat$S2,dat$medTC,log="y")

#make correlation plots
plot(dat$S1[which(dat$cor<0)],dat$E1[which(dat$cor<0)],pch='-',cex=1.5, col=ifelse(dat$CoexistRank[which(dat$cor<0)]==0,"black",ifelse(dat$CoexistRank[which(dat$cor<0)]==1,"yellow","green")),ylim=c(0.5,1),xlab="Strength of stabilization",ylab="Fitness equivalence")
identify(dat$S1[which(dat$cor<0)],dat$E1[which(dat$cor<0)],labels=dat$ID[which(dat$cor<0)])
plot(dat$S1[which(dat$cor>=0)],dat$E1[which(dat$cor>=0)],pch='+',cex=1.5, col=ifelse(dat$CoexistRank[which(dat$cor>=0)]==0,"black",ifelse(dat$CoexistRank[which(dat$cor>=0)]==1,"yellow","green")),ylim=c(0.5,1),xlab="Strength of stabilization",ylab="Fitness equivalence")
identify(dat$S1[which(dat$cor>=0)],dat$E1[which(dat$cor>=0)],labels=dat$ID[which(dat$cor>=0)])

#make boxplot of coexistence vs strength of stabil by rank
boxplot(S1~CoexistRank,data=dat, names=c("1","1","1"), col=c("white","lightgrey","darkgrey"), xlab="Species", ylab="Strength of stabilization",at=c(1,4,7),xlim=c(0,9),ylim=c(1,7.8),range=0)
boxplot(S2~CoexistRank,data=dat, names=c("2","2","2"), col=c("white","lightgrey","darkgrey"), add=T,at=c(2,5,8),range=0)
abline(v=3,lwd=2)
abline(v=6,lwd=2)
text(1.2,7.7,labels="Persistence",cex=1)
text(1.4,7.2,labels="<100",cex=1)
text(4.5,7.7,labels="100<t<1000",cex=1)
text(7.4,7.7,labels=">1000",cex=1)

#plot of persistence vs asymmetry of stabilization
plot(dat$Asy,dat$medTC,log="y",xlab="difference in SOS", ylab="Median persistance time")

dat$Asycat=NA
dat$Asycat[which(dat$Asy>0)]=1
dat$Asycat[which(dat$Asy==0)]=0
dat$Asycat[which(dat$Asy<0)]=-1
dat$Asycat[which(dat$Asy<=-2)]=-2
dat$Asycat[which(dat$Asy>=2)]=2
dat$Asycat[which(dat$Asy<=-3)]=-3
dat$Asycat[which(dat$Asy>=3)]=3

boxplot(medTC~Asycat,data=dat,xlab="Difference in stabilization",ylab="Median persistence time",outline=F)

#log-linear model of the effect of stabilization on persistence
summary(lm(log(medTC)~S1*E1+cor,data=dat))

coexist2=dat[which(dat$CoexistRank==2),]
coexist1=dat[which(dat$CoexistRank==1),]
extinct=dat[which(dat$CoexistRank==0),]
null=t.test(dat$cor)
coex2=t.test(coexist2$cor)
coex1=t.test(coexist1$cor)
ext=t.test(extinct$cor)
library(Hmisc)
errbar(c("<100","all","100<t<1000",">1000"),c(ext$estimate,null$estimate,coex1$estimate,coex2$estimate),c(ext$conf.int[2],null$conf.int[2],coex1$conf.int[2],coex2$conf.int[2]),c(ext$conf.int[1],null$conf.int[1],coex1$conf.int[1],coex2$conf.int[1]),ylab="Relationship between SOS and abundance")



#effect of l1, l2
plot(dat$l1,dat$l2,col=ifelse(dat$CoexistRank==0,"black",ifelse(dat$CoexistRank==1,"yellow","green")))
layout(matrix(c(1,2,3),1,3))
plot(dat$l1[which(dat$CoexistRank==2)]+runif(length(which(dat$CoexistRank==2))),dat$l2[which(dat$CoexistRank==2)]+runif(length(which(dat$CoexistRank==2))),xlim=c(10,21),ylim=c(10,21),col="black",pch=19,xlab="l1",ylab="l2",main="t>1000")
plot(dat$l1[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1))),dat$l2[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1))),col="grey",pch=19,xlab="l1",ylab="l2",main="100<t<1000")
plot(dat$l1[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0))),dat$l2[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0))),col="grey",xlab="l1",ylab="l2",main="t<100")



