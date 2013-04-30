#Stochastic Lotka-Volterra logistic growth model

# Parameters--------------------------------

r=c(1.5,1.5,1.5)
K=c(100,80,50)
TotTime=800
initialNsp1 = 5
initialNsp2 = 5
initialNsp3 = 5

a11=1
a12=.3
a13=.1
a22=1
a21=.3
a23=.3
a33=1
a31=.1
a32=.3

zvar1=.02
zvar2=.02
zvar3=.02
zcov=-0.01
zmeans = c(0,0,0)

#-------------------------------------------

library(mvtnorm)

T=1:TotTime

updateN = function(r_arg, K_arg, Nself, N1, N2, a_intra, a1, a2, z_arg){
    newN = Nself + Nself * (r_arg * ((K_arg - a_intra*
    Nself - a1 * N1 - a2*N2) / K_arg) + z_arg)
    if(newN < 0) newN = 0
    return(newN)
}

N = matrix(data=NA,nrow=TotTime,ncol=3)
N[1,1] = initialNsp1
N[1,2] = initialNsp2
N[1,3] = initialNsp3

#generate the sigma matrix for generation of random z's
sigma = matrix(data=c(zvar1,zcov,zcov,zcov,zvar2,zcov,zcov,zcov,zvar3),nrow=3,ncol=3)
z = rmvnorm(TotTime,zmeans,sigma)


# main loop ---------------------------------

for(i in 2:TotTime) {
    N[i,1] = floor(updateN(r[1],K[1],N[i-1,1],N[i-1,2],N[i-1,3],a11,a12,a13,z[i,1]))
    N[i,2] = floor(updateN(r[2],K[2],N[i-1,2],N[i-1,1],N[i-1,3],a22,a21,a23,z[i,2]))
    N[i,3] = floor(updateN(r[3],K[3],N[i-1,3],N[i-1,1],N[i-1,2],a33,a31,a32,z[i,3]))
}

plot(N[,1],type="o",xlab="Time", col="steelblue", ylab="N",ylim=c(0,180))
par(new=TRUE)
plot(N[,2],type="o",xlab="Time", col="firebrick", ylab="",ylim=c(0,180))
par(new=TRUE)
plot(N[,3],type="o",xlab="Time", col="lightgreen", ylab="",ylim=c(0,180))

#calculate observed growth rates:

#log lambda function
lambda=function(file) {
if(file[i+1,x]==0 | file[i,x]==0){
  r=NA
}else{
  r=log(file[i+1,x]/file[i,x])
}
return(r)
}

#calc frequency
C=apply(N, 1, sum)
freq=N/C

#calc lambda
rates=matrix(data=NA,nrow=length(N[,1])-1,ncol=length(N[1,]))
for(i in 1:(length(N[,1])-1)) {
  for(x in 1:length(N[1,])) {
    rates[i,x]=lambda(N)
  }
}


plot(freq[1:length(freq[,1])-1,1],rates[,1],xlim=c(0,1),ylim=c(-3,3),col="steelblue",xlab="Frequency",ylab="log r")
points(freq[1:length(freq[,2])-1,2],rates[,2],pch=2,col="firebrick")
points(freq[1:length(freq[,3])-1,3],rates[,3],pch=3,col="lightgreen")
abline(h=0,col="grey")

abline(lm(rates[,1]~freq[1:length(freq[,1])-1,1]),col="steelblue")
abline(lm(rates[,2]~freq[1:length(freq[,2])-1,2]),col="firebrick")
abline(lm(rates[,3]~freq[1:length(freq[,3])-1,3]),col="lightgreen")

