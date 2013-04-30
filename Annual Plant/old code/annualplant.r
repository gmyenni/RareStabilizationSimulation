#Adler et al annual plant growth model

# Parameters--------------------------------

TotTime=200
#l1=20
#l2=15
#l3=10

initialNsp1 = 5
initialNsp2 = 5
initialNsp3 = 5

a11=1
#a12=.7
#a13=.2
a22=1
#a21=.7
#a23=.2
a33=1
#a31=.2
#a32=.2


#-------------------------------------------
#generate random lambdas
r1=rpois(TotTime,l1)
r2=rpois(TotTime,l2)
r3=rpois(TotTime,l3)

T=1:TotTime

updateN = function(r_arg, Nself, N1, N2, a_intra, a1, a2){
    newN = r_arg * Nself /(1 + a_intra * Nself + a1 * N1 + a2 * N2)
    if(newN < 0) newN = 0
    return(newN)
}

N = matrix(data=NA,nrow=TotTime,ncol=3)
N[1,1] = initialNsp1
N[1,2] = initialNsp2
N[1,3] = initialNsp3



# main loop ---------------------------------

for(t in 2:TotTime) {
    N[t,1] = floor(updateN(r1[t-1],N[t-1,1],N[t-1,2],N[t-1,3],a11,a12,a13))
    N[t,2] = floor(updateN(r2[t-1],N[t-1,2],N[t-1,1],N[t-1,3],a22,a21,a23))
    N[t,3] = floor(updateN(r3[t-1],N[t-1,3],N[t-1,1],N[t-1,2],a33,a31,a32))
}

plot(N[,1],type="o",xlab="Time", col="steelblue", ylab="N",ylim=c(0,30))
par(new=TRUE)
plot(N[,2],type="o",xlab="Time", col="firebrick", ylab="",ylim=c(0,30))
par(new=TRUE)
plot(N[,3],type="o",xlab="Time", col="lightgreen", ylab="",ylim=c(0,30))

##calculate observed growth rates:
#
##log lambda function
#lambda=function(file) {
#if(file[i+1,x]==0 | file[i,x]==0){
#  r=NA
#}else{
#  r=log(file[i+1,x]/file[i,x])
#}
#return(r)
#}
#
#calc frequency
#C=apply(N, 1, sum)
#freq=N/C
#
#calc lambda
#rates=matrix(data=NA,nrow=length(N[,1])-1,ncol=length(N[1,]))
#for(i in 1:(length(N[,1])-1)) {
#  for(x in 1:length(N[1,])) {
#    rates[i,x]=lambda(N)
#  }
#}
#
#
#plot(freq[1:length(freq[,1])-1,1],rates[,1],xlim=c(0,1),ylim=c(-3,3),col="steelblue",xlab="Frequency",ylab="log r")
#points(freq[1:length(freq[,2])-1,2],rates[,2],pch=2,col="firebrick")
#points(freq[1:length(freq[,3])-1,3],rates[,3],pch=3,col="lightgreen")
#abline(h=0,col="grey")
#
#abline(lm(rates[,1]~freq[1:length(freq[,1])-1,1]),col="steelblue")
#abline(lm(rates[,2]~freq[1:length(freq[,2])-1,2]),col="firebrick")
#abline(lm(rates[,3]~freq[1:length(freq[,3])-1,3]),col="lightgreen")
#
#cor.test(freq[1:length(freq[,1])-1,1],rates[,1])
#cor.test(freq[1:length(freq[,2])-1,2],rates[,2])
#cor.test(freq[1:length(freq[,3])-1,3],rates[,3])