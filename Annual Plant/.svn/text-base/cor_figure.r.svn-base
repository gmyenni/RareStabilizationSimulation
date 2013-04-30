dat_det=read.csv("annplant_2spp_det_data.txt")
dat_stoch=read.csv("annplant_2spp_stoch_data.txt")

dat_det=dat_det[which(dat_det$Rank1==2),]
dat_stoch=dat_stoch[which(dat_stoch$Rank1==2),]

#calculate correlation btw abund and stabilization
corrs=rep(NA,length(dat_det[,1]))
for(i in 1:length(dat_det[,1])) {
corrs[i]=cov(x=c(dat_det[i,6],dat_det[i,8]),y=c(dat_det[i,16],dat_det[i,18]))  }
dat_det$cor=corrs

#calculate correlation btw abund and stabilization
corrs=rep(NA,length(dat_stoch[,1]))
for(i in 1:length(dat_stoch[,1])) {
corrs[i]=cov(x=c(dat_stoch[i,6],dat_stoch[i,8]),y=c(dat_stoch[i,21],dat_stoch[i,23]))  }
dat_stoch$cor=corrs


#make correlation plots
layout(matrix(c(1,3,2,4),2,2))

plot(dat_det$S1[which(dat_det$cor<0)],dat_det$E1[which(dat_det$cor<0)],cex=1.5, pch=ifelse(dat_det$coexist[which(dat_det$cor<0)]==500,19,1),col=ifelse(dat_det$coexist[which(dat_det$cor<0)]==500,"black","lightgrey"),ylim=c(0.5,1),xlab=NA,ylab="Fitness equivalence",main="Negative Correlation")

plot(dat_det$S1[which(dat_det$cor>=0)],dat_det$E1[which(dat_det$cor>=0)],pch=ifelse(dat_det$coexist[which(dat_det$cor>=0)]==500,19,1),cex=1.5, col=ifelse(dat_det$coexist[which(dat_det$cor>=0)]==500,"black","lightgrey"),ylim=c(0.5,1),xlab=NA,ylab=NA,main="Positive Correlation")


plot(dat_stoch$S1[which(dat_stoch$cor<0)],dat_stoch$E1[which(dat_stoch$cor<0)],pch=ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor<0)]==0,1,ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor<0)]==1,19,19)),cex=1.5, col=ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor<0)]==0,"lightgrey",ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor<0)]==1,"darkgrey","black")),ylim=c(0.5,1),xlab="Strength of stabilization",ylab="Fitness equivalence")

plot(dat_stoch$S1[which(dat_stoch$cor>=0)],dat_stoch$E1[which(dat_stoch$cor>=0)],pch=ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor>=0)]==0,1,ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor>=0)]==1,19,19)),cex=1.5, col=ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor>=0)]==0,"lightgrey",ifelse(dat_stoch$CoexistRank[which(dat_stoch$cor>=0)]==1,"darkgrey","black")),ylim=c(0.5,1),xlab="Strength of stabilization",ylab=NA)
