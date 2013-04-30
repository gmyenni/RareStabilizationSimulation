#Adler et al annual plant growth model(deterministic)

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
    N[t,1] = floor(updateN(l1,N[t-1,1],N[t-1,2],N[t-1,3],a11,a12,a13))
    N[t,2] = floor(updateN(l2,N[t-1,2],N[t-1,1],N[t-1,3],a22,a21,a23))
    N[t,3] = floor(updateN(l3,N[t-1,3],N[t-1,1],N[t-1,2],a33,a31,a32))
}

#plot(N[,1],type="o",xlab="Time", col="steelblue", ylab="N",ylim=c(0,30))
#par(new=TRUE)
#plot(N[,2],type="o",xlab="Time", col="firebrick", ylab="",ylim=c(0,30))
#par(new=TRUE)
#plot(N[,3],type="o",xlab="Time", col="lightgreen", ylab="",ylim=c(0,30))
