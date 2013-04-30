#Stochastic logistic growth model, 4 spp, env stochasticity

# Parameters--------------------------------

TotTime=200
initialNsp1 = 20
initialNsp2 = 20
initialNsp3 = 20
initialNsp4 = 20
totspp=4
K=100            #community carrying capacity

#frequency dependence
sp1=-1                  #species 1 slope
sp2=-1
sp3=-2
sp4=-4
int1=.8                 #species 1 intercept
int2=.5
int3=.3
int4=.15

#Demographic stochasticity parameters
zvar1=.03
zvar2=.03
zvar3=.03
zvar4=.03
zcov=0
zmeans = c(0,0,0,0)

#Environmental stochasticity parameters
evar1=.04
evar2=.04
evar3=.04
evar4=.04
ecov=-.05
emeans = c(1,1,1,1)
#-------------------------------------------

library(mvtnorm)

T=1:TotTime

#frequency dependence model
rdep=function(slope,int,N,Ctot){
    if(Ctot==N) r=1 else r=exp(slope*(N/Ctot-int))
    return(r)
}

#population model
updateN = function(r_arg, z_arg, e_arg, N, Ctot){
    newN = (r_arg + z_arg)* N * (K*e_arg)/Ctot
    if(newN < 0) newN = 0
    return(newN)
}

N = matrix(data=NA,nrow=TotTime,ncol=totspp)
N[1:2,1] = initialNsp1
N[1:2,2] = initialNsp2
N[1:2,3] = initialNsp3
N[1:2,4] = initialNsp4

#generate the sigma matrices for generation of random z's and e's
sigma = matrix(data=c(zvar1,zcov,zcov,zcov,zcov,zvar2,zcov,zcov,zcov,zcov,zvar3,zcov,zcov,zcov,zcov,zvar4),nrow=4,ncol=4)
z = rmvnorm(TotTime,zmeans,sigma)

sigma_env = matrix(data=c(evar1,ecov,ecov,ecov,ecov,evar2,ecov,ecov,ecov,ecov,evar3,ecov,ecov,ecov,ecov,evar4),nrow=4,ncol=4)
e = rmvnorm(TotTime,emeans,sigma_env)

# main loop ---------------------------------

for(i in 3:TotTime) {
    N[i,1] = floor(updateN(rdep(sp1,int1,N[i-2,1],sum(N[i-2,])),z[i,1],e[i,1],N[i-1,1],sum(N[i-1,])))
    N[i,2] = floor(updateN(rdep(sp2,int2,N[i-2,2],sum(N[i-2,])),z[i,2],e[i,2],N[i-1,2],sum(N[i-1,])))
    N[i,3] = floor(updateN(rdep(sp3,int3,N[i-2,3],sum(N[i-2,])),z[i,3],e[i,3],N[i-1,3],sum(N[i-1,])))
    N[i,4] = floor(updateN(rdep(sp4,int4,N[i-2,4],sum(N[i-2,])),z[i,4],e[i,4],N[i-1,4],sum(N[i-1,])))
}


plot(N[,1],type="o",xlab="", ylab="", xaxt='n', yaxt='n', col="steelblue", ylim=c(0,max(N)))
points(N[,2],type="o", col="firebrick")
points(N[,3],type="o", col="lightgreen")
points(N[,4],type="o", col="orange")
axis(side=3,xlab="Time")
axis(side=4,ylab="N")
mtext("N", side=4)
mtext("Time", side=3)
par(new=TRUE)
plot((0:K)/K,log(rdep(sp1,int1,0:K,K)),ylab="log r",xlab="Frequency",pch=19,ylim=c(-3,3), col="steelblue")
points((0:K)/K,log(rdep(sp2,int2,0:K,K)),pch=19, col="firebrick")
points((0:K)/K,log(rdep(sp3,int3,0:K,K)),pch=19, col="lightgreen")
points((0:K)/K,log(rdep(sp4,int4,0:K,K)),pch=19, col="orange")
abline(h=0,col="grey")


#calculate observed growth rates:

#log lambda function
#lambda=function(file) {
#if(file[i+1,x]==0 | file[i,x]==0){
#  r=NA
#}else{
#  r=log(file[i+1,x]/file[i,x])
#}
#return(r)
#}

#calc frequency
#C=apply(N, 1, sum)
#freq=N/C

#calc lambda
#rates=matrix(data=NA,nrow=length(N[,1])-1,ncol=length(N[1,]))
#for(i in 1:(length(N[,1])-1)) {
#  for(x in 1:length(N[1,])) {
#    rates[i,x]=lambda(N)
#  }
#}


#plot(freq[1:length(freq[,1])-1,1],rates[,1],xlim=c(0,1),ylim=c(-3,3),col="steelblue",xlab="Frequency",ylab="log r")
#points(freq[1:length(freq[,2])-1,2],rates[,2],pch=2,col="firebrick")
#points(freq[1:length(freq[,3])-1,3],rates[,3],pch=3,col="lightgreen")
#points(freq[1:length(freq[,4])-1,4],rates[,4],pch=4,col="orange")
#abline(h=0,col="grey")

#abline(lm(rates[,1]~freq[1:length(freq[,1])-1,1]),col="steelblue")
#abline(lm(rates[,2]~freq[1:length(freq[,2])-1,2]),col="firebrick")
#abline(lm(rates[,3]~freq[1:length(freq[,3])-1,3]),col="lightgreen")
#abline(lm(rates[,4]~freq[1:length(freq[,4])-1,4]),col="orange")
