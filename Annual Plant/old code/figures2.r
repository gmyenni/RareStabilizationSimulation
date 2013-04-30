dat1=read.csv("annplant_stat2")
dat1=dat1[2:length(dat1[,1]),]
dat2=read.csv("annplant_stat3")
dat2=dat2[2:length(dat2[,1]),]

data1=rbind(dat1,dat2)
data2=data1[,2:length(data1[1,])]

#calculate strengths of stabilization
#S=gm(l1,l2)*(a12 a21 - a11 a22)/(a12 a21 + a11 a32 + a22 a31 - a11 a22 - a12 a31 - a21 a32)+(a12 a31 - a11 a32)l2+(a21 a32 - a22 a31)l1   (for species 3)

S1=function(S){
    l1=S[1]
    l2=S[2]
    l3=S[3]
    a12=S[4]
    a13=S[5]
    a21=S[6]
    a23=S[7]
    a31=S[8]
    a32=S[9]
    s=mean(l2,l3)*(a32*a23 - 1)/(a32*a23 + a12 + a13 - 1 - a32*a13 - a23*a12)+(a32*a13 - a12)*l2+(a23*a12 - a13)*l3
    return(s)
}

S2=function(S){
    l1=S[1]
    l2=S[2]
    l3=S[3]
    a12=S[4]
    a13=S[5]
    a21=S[6]
    a23=S[7]
    a31=S[8]
    a32=S[9]
    s=mean(l1,l3)*(a13*a31 - 1)/(a13*a31 + a23 + a21 - 1 - a13*a21 - a31*a23)+(a13*a21 - a23)*l3+(a31*a23 - a21)*l1
    return(s)
}

S3=function(S){
    l1=S[1]
    l2=S[2]
    l3=S[3]
    a12=S[4]
    a13=S[5]
    a21=S[6]
    a23=S[7]
    a31=S[8]
    a32=S[9]
    s=mean(l1,l2)*(a12*a21 - 1)/(a12*a21 + a32 + a31 - 1 - a12*a31 - a21*a32)+(a12*a31 - a32)*l2+(a21*a32 - a31)*l1
    return(s)
}
    
#calculate fitness inequivalence
#FI=l3/gm(l1,l2)     (for species 3)

FI1=function(L){
  l1=L[1]
  l2=L[2]
  l3=L[3]
  fi=l1/mean(l2,l3)
  return(fi)
}

FI2=function(L){
  l1=L[1]
  l2=L[2]
  l3=L[3]
  fi=l2/mean(l1,l3)
  return(fi)
}

FI3=function(L){
  l1=L[1]
  l2=L[2]
  l3=L[3]
  fi=l3/mean(l1,l2)
  return(fi)
}

#fitness - stabilization matrix
FS=cbind(apply(data2[,1:3],1,FI1),apply(data2[,1:3],1,FI2),apply(data2[,1:3],1,FI3),apply(data2[,1:9],1,S1),apply(data2[,1:9],1,S2),apply(data2[,1:9],1,S3))

#calc evenness index of fitness equivalence & stabilization
library(vegan)
FS2=cbind(FS,diversity(FS[,1:3], index = "shannon", MARGIN = 1),diversity(FS[,4:6], index = "shannon", MARGIN = 1))

#calculate correlation btw abund and stabilization
corrs=rep(NA,length(data2[,1]))


for(i in 1:length(data2[,1])) {

corrs[i]=cov(x=c(data2[i,10],data2[i,11],data2[i,12]),y=c(FS[i,4],FS[i,5],FS[i,6]))  }


s=which(apply(ifelse(data2[,13:15]>20,1,0),1,sum)!=3)
t=which(apply(ifelse(data2[,13:15]>20,1,0),1,sum)==3)

plot(1/FS2[,8],FS2[,7],pch=ifelse(corrs>0,'+','-'),xlim=c(0.9,1.05),col=ifelse(apply(ifelse(data2[,13:15]>20,1,0),1,sum)==3,'green',rgb(0, 0, 0,0.1)),cex=1.5,xlab="Strength of stabilization",ylab="Fitness equivalence")                                           

plot(1/FS2[s,8],FS2[s,7],pch=ifelse(corrs[s]>0,'+','-'),xlim=c(0.9,1.05),col='black',xlab="Strength of stabilization",ylab="Fitness equivalence")                                           

points(1/FS2[t,8],FS2[t,7],pch=ifelse(corrs[t]>0,'+','-'),col='green',cex=1.5)  