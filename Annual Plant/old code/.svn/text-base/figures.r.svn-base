dat=read.csv("annplant_stat")
dat=dat[2:length(dat[,1]),1:22]
dat1=read.csv("annplant_stat10")
dat1=dat1[2:length(dat1[,1]),]
dat2=read.csv("annplant_stat11")
dat2=dat2[2:length(dat2[,1]),]
dat3=read.csv("annplant_stat12")
dat3=dat3[2:length(dat3[,1]),]
dat4=read.csv("annplant_stat102")
dat4=dat4[2:length(dat4[,1]),]
dat5=read.csv("annplant_stat112")
dat5=dat5[2:length(dat5[,1]),]
dat6=read.csv("annplant_stat122")
dat6=dat6[2:length(dat6[,1]),]
dat7=read.csv("annplant_stat103")
dat7=dat7[2:length(dat7[,1]),]
dat8=read.csv("annplant_stat113")
dat8=dat8[2:length(dat8[,1]),]
dat9=read.csv("annplant_stat123")
dat9=dat9[2:length(dat9[,1]),]
dat10=read.csv("annplant_stat132")
dat10=dat10[2:length(dat10[,1]),]
dat11=read.csv("annplant_stat133")
dat11=dat11[2:length(dat11[,1]),]
dat12=read.csv("annplant_stat143")
dat12=dat12[2:length(dat12[,1]),]
dat13=read.csv("annplant_stat105")
dat13=dat13[2:length(dat13[,1]),]
dat14=read.csv("annplant_stat106")
dat14=dat14[2:length(dat14[,1]),]
dat15=read.csv("annplant_stat107")
dat15=dat15[2:length(dat15[,1]),]
dat16=read.csv("annplant_stat108")
dat16=dat16[2:length(dat16[,1]),]
dat17=read.csv("annplant_stat109")
dat17=dat17[2:length(dat17[,1]),]
dat18=read.csv("annplant_stat119")
dat18=dat18[2:length(dat18[,1]),]

data1=rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11,dat12,dat13,dat14,dat15,dat16,dat17,dat18)
data2=data1[,2:length(data1[1,])]
data2[,16]=ifelse(is.na(data2[,16])==T,200,data2[,16])
data2[,17]=ifelse(is.na(data2[,17])==T,200,data2[,17])
data2[,18]=ifelse(is.na(data2[,18])==T,200,data2[,18])

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

#calc richness index of fitness equivalence and stabilization
library(vegan)
FS2=cbind(FS,diversity(FS[,1:3], index = "shannon", MARGIN = 1),diversity(FS[,4:6], index = "shannon", MARGIN = 1))

#calculate correlation btw abund and stabilization
corrs=rep(NA,length(data2[,1]))


for(i in 1:length(data2[,1])) {

corrs[i]=cov(x=c(data2[i,10],data2[i,11],data2[i,12]),y=c(FS[i,4],FS[i,5],FS[i,6]))  }


s=which(apply(ifelse(data2[,16:18]>20,1,0),1,sum)!=3)
t=which(apply(ifelse(data2[,16:18]>20,1,0),1,sum)==3)

plot(1/FS2[,8],FS2[,7],pch=ifelse(corrs>0,'+','-'),xlim=c(0.9,1.05),col=ifelse(apply(ifelse(data2[,19:21]>5,1,0),1,sum)==3,'green',rgb(0, 0, 0,0.1)),cex=1.5,xlab="Strength of stabilization",ylab="Fitness equivalence")                                           

plot(1/FS2[s,8],FS2[s,7],pch=ifelse(corrs[s]>0,'+','-'),xlim=c(0.9,1.05),col='black',xlab="Strength of stabilization",ylab="Fitness equivalence")                                           

points(1/FS2[t,8],FS2[t,7],pch=ifelse(corrs[t]>0,'+','-'),col='green',cex=1.5)                                           
