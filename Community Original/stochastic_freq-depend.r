#Stochastic neg-freq depend logistic growth model, 3 species

# Parameters--------------------------------

TotTime=200
initialNsp1 = 30
initialNsp2 = 30
initialNsp3 = 30
totspp=3
K=100             #community carrying capacity

#frequency dependence
#sp1=-1                  #species 1 slope
#sp2=-1
#sp3=-1
#int1=.6                 #species 1 intercept
#int2=.4
#int3=.1

zvar1=.01
zvar2=.01
zvar3=.01
zcov=-.01
zmeans = c(0,0,0)

#-------------------------------------------

library(mvtnorm)

T=1:TotTime

#frequency dependence model
rdep=function(slope,int,N,Ctot){
    if(Ctot==N) r=1 else r=exp(slope*(N/Ctot-int))
    return(r)
}

#population model
updateN = function(r_arg, z_arg, N, Ctot){
    newN = (r_arg + z_arg)* N * K/Ctot
    if(newN < 0) newN = 0
    return(newN)
}

N = matrix(data=NA,nrow=TotTime,ncol=totspp)
N[1:2,1] = initialNsp1
N[1:2,2] = initialNsp2
N[1:2,3] = initialNsp3

#generate the sigma matrix for generation of random z's
sigma = matrix(data=c(zvar1,zcov,zcov,zcov,zvar2,zcov,zcov,zcov,zvar3),nrow=3,ncol=3)
z = rmvnorm(TotTime,zmeans,sigma)


# main loop ---------------------------------

for(i in 3:TotTime) {
    N[i,1] = floor(updateN(rdep(sp1,int1,N[i-1,1],sum(N[i-1,])),z[i,1],N[i-1,1],sum(N[i-1,])))
    N[i,2] = floor(updateN(rdep(sp2,int2,N[i-1,2],sum(N[i-1,])),z[i,2],N[i-1,2],sum(N[i-1,])))
    N[i,3] = floor(updateN(rdep(sp3,int3,N[i-1,3],sum(N[i-1,])),z[i,3],N[i-1,3],sum(N[i-1,])))
  
}


plot(N[,1],type="o",xlab="", ylab="", xaxt='n', yaxt='n', col="steelblue", ylim=c(0,max(N)))
points(N[,2],type="o", col="firebrick")
points(N[,3],type="o", col="lightgreen")
axis(side=3,xlab="Time")
axis(side=4,ylab="N")
mtext("N", side=4)
mtext("Time", side=3)
par(new=TRUE)
plot((0:K)/K,log(rdep(sp1,int1,0:K,K)),ylab="log r",xlab="Frequency",pch=19,ylim=c(-3,3), col="steelblue")
points((0:K)/K,log(rdep(sp2,int2,0:K,K)),pch=19, col="firebrick")
points((0:K)/K,log(rdep(sp3,int3,0:K,K)),pch=19, col="lightgreen")
abline(h=0,col="grey")