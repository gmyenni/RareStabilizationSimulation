#calculate relationship between abundance and stabilization for all simulation results
#species 1 is the inferior species
dat=read.csv("annplant_2spp_det_rare.txt")

#calculate correlation btw abund and stabilization
#corrs=rep(NA,length(dat[,1]))
#for(i in 1:length(dat[,1])) {
#corrs[i]=cov(x=c(dat[i,8],dat[i,9]),y=c(dat[i,11],dat[i,13]))  }
#dat$cor=corrs
#write.csv(dat,"annplant_2spp_det.csv")

#make correlation plots
#plot(dat$S1,dat$E1,pch=ifelse(dat$cor>0,'+','-'),col=ifelse(dat$coexist==50,"green","black"),xlab="Strength of stabilization",ylab="Fitness equivalence")

plot(dat$S1[which(dat$cor<0)],dat$E1[which(dat$cor<0)],pch='-',cex=1.5, col=ifelse(dat$CoexistRank[which(dat$cor<0)]==1,"green","black"),xlab="Strength of stabilization",ylab="Fitness equivalence")

plot(dat$S1[which(dat$cor>=0)],dat$E1[which(dat$cor>=0)],pch='+',cex=1.5, col=ifelse(dat$CoexistRank[which(dat$cor>=0)]==1,"green","black"),xlab="Strength of stabilization",ylab="Fitness equivalence")

#effect on coexistence
summary(glm(CoexistRank~S1*E1+cor,data=dat,family=binomial(link="logit")))
t.test(cor~CoexistRank,data=dat)
coexist=dat[which(dat$CoexistRank==1),]
extinct=dat[which(dat$CoexistRank==0),]
null=t.test(dat$cor)
coex=t.test(coexist$cor)
ext=t.test(extinct$cor)
library(Hmisc)
errbar(c("extinct","all","coexist"),c(ext$estimate,null$estimate,coex$estimate),c(ext$conf.int[2],null$conf.int[2],coex$conf.int[2]),c(ext$conf.int[1],null$conf.int[1],coex$conf.int[1]),ylab="Relationship between SOS and abundance")

#make boxplot of coexistence vs strength of stabil by rank
boxplot(S1~CoexistRank,data=dat, names=c("1","1"),col=c("white","grey"),xlab="Species", ylab="Strength of Stabilization",at=c(1,4),xlim=c(0,6),ylim=c(1,13),range=0)
boxplot(S2~CoexistRank,data=dat, names=c("2","2"),add=T,at=c(2,5),col=c("white","grey"),range=0)
abline(v=3,lwd=2)
text(1.4,7.7,labels="Competitive",cex=1)
text(1.4,7.2,labels="exclusion",cex=1)
text(4.63,7.7,labels="Coexistence",cex=1)

#effect of l1, l2
layout(matrix(c(1,2),1,2))
plot(dat$l1[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1))),dat$l2[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1))),xlim=c(10,21),ylim=c(10,21),col="black",pch=19,xlab="l1",ylab="l2",main="Coexistence")
plot(dat$l1[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0))),dat$l2[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0))),col="grey",xlab="l1",ylab="l2",main="Comp. Exclusion")

#effect of a11, a22
layout(matrix(c(1,2),1,2))
plot(dat$a11[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a22[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,col="black",pch=19,xlab="a11",ylab="a22",main="Coexistence")
plot(dat$a11[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a22[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,col="grey",xlab="a11",ylab="a22",main="Comp. Exclusion")

#effect of a12, a21
layout(matrix(c(1,2),1,2))
plot(dat$a12[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a21[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,xlim=c(0,1.1),ylim=c(0,1.1),col="black",pch=19,xlab="a12",ylab="a21",main="Coexistence")
plot(dat$a12[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a21[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,xlim=c(0,1.1),ylim=c(0,1.1),col="grey",xlab="a12",ylab="a21",main="Comp. Exclusion")

#effect of a12, a22
layout(matrix(c(1,2),1,2))
plot(dat$a12[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a22[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,col="black",pch=19,xlab="a12",ylab="a22",main="Coexistence")
plot(dat$a12[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a22[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,col="grey",xlab="a12",ylab="a22",main="Comp. Exclusion")

#effect of a21, a11
layout(matrix(c(1,2),1,2))
plot(dat$a21[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a11[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,col="black",pch=19,xlab="a21",ylab="a11",main="Coexistence")
plot(dat$a21[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a11[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,col="grey",xlab="a21",ylab="a11",main="Comp. Exclusion")

#effect of E1, a11
layout(matrix(c(1,2),1,2))
plot(dat$E1[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a11[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,col="black",pch=19,xlab="E1",ylab="a11",main="Coexistence")
plot(dat$E1[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a11[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,col="grey",xlab="E1",ylab="a11",main="Comp. Exclusion")

#effect of E2, a22
layout(matrix(c(1,2),1,2))
plot(dat$E2[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a22[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,col="black",pch=19,xlab="E2",ylab="a22",main="Coexistence")
plot(dat$E2[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a22[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,col="grey",xlab="E2",ylab="a22",main="Comp. Exclusion")

#effect of E1, a22
layout(matrix(c(1,2),1,2))
plot(dat$E1[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a22[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,col="black",pch=19,xlab="E1",ylab="a22",main="Coexistence")
plot(dat$E1[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a22[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,col="grey",xlab="E1",ylab="a22",main="Comp. Exclusion")

#effect of E2, a11
layout(matrix(c(1,2),1,2))
plot(dat$E2[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,dat$a11[which(dat$CoexistRank==1)]+runif(length(which(dat$CoexistRank==1)))/10,col="black",pch=19,xlab="E2",ylab="a11",main="Coexistence")
plot(dat$E2[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,dat$a11[which(dat$CoexistRank==0)]+runif(length(which(dat$CoexistRank==0)))/10,col="grey",xlab="E2",ylab="a11",main="Comp. Exclusion")
