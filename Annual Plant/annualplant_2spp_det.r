print(date())
# parameters -----------------------------------------
  # name of output file:
  outfile =  "annplant_2spp_det.csv"

  #specify combinations of freq-dep parameters
  l1_v=10:20
  l2_v=10:20

  a11_v=seq(0.1,1,by=0.1)
  a12_v=seq(0.1,1,by=0.1)
  a21_v=seq(0.1,1,by=0.1)
  a22_v=seq(0.1,1,by=0.1)

  simul = expand.grid("l1"=l1_v,"l2"=l2_v,"a11"=a11_v,"a12"=a12_v,"a21"=a21_v,"a22"=a22_v)


#2-species annual plant model analytical solution----------------------------------------------
analyN = function(r1, r2, a1, a12, a21, a2){
    N1 = (r1-1-(a12/a2)*(r2-1))/(a1-a21*a12/a2)
    N2 = (r2-1-(a21/a1)*(r1-1))/(a2-a21*a12/a1)
    return(c(N1,N2))
}

#simulation output
simul$N1 = NA
simul$N2 = NA
simul$coexist = NA

#simulation loop-------------------------------------
for(iSim in 1:dim(simul)[1]) {

    l1=simul$l1[iSim]           #get parameter combo
    l2=simul$l2[iSim]
    a11=simul$a11[iSim]
    a12=simul$a12[iSim]
    a21=simul$a21[iSim]
    a22=simul$a22[iSim]
    
#simulation stats
simul$N1[iSim] = analyN(l1,l2,a11,a12,a21,a22)[1]
simul$N2[iSim] = analyN(l1,l2,a11,a12,a21,a22)[2]
simul$coexist[iSim] = ifelse(sum((simul[iSim,7:8]<=0)>=1),0,1)

#progress update
  tmp = paste("Parameter combination",iSim,"of",dim(simul)[1],sep=" ")
  print(tmp)
  print(date())
  flush.console()
} #next parm combo (iSim)


write.table(simul,outfile, sep=",",row.names=FALSE)


# if only one run performed, make figures
if(dim(simul)[1]==1){
  par(mfrow=c(1,2),tcl=-0.2)
  # plot density time series
  matplot(N,type="l",xlab="Time",ylab="N",col=c("steelblue","firebrick"))
  # plot frequency dependence in growth
  Nfreq = N/rowSums(N)  # calculate frequencies
  Nfreq[Nfreq==1]=NA # remove values after one species goes exinct
  growth = log(N[2:NROW(N),])-log(N[1:(NROW(N)-1),])
  growth[growth==-Inf]=NA
  myLims = c(min(growth,na.rm=T)-0.05,max(growth,na.rm=T)+0.05)
  plot(Nfreq[1:NROW(Nfreq)-1,1],growth[,1],xlab="Frequency",ylab="Growth rate",
    xlim=c(0,1),ylim=myLims,col="steelblue")
  abline(lm(growth[,1]~ Nfreq[1:NROW(Nfreq)-1,1] ),col="steelblue")
  par(new=T)
  plot(Nfreq[1:NROW(Nfreq)-1,2],growth[,2],xlab="",ylab="",xaxt="n",yaxt="n",
    xlim=c(0,1),ylim=myLims,col="firebrick")
  abline(lm(growth[,2]~ Nfreq[1:NROW(Nfreq)-1,2] ),col="firebrick")
  abline(0,0,lty="dotted")
  legend("topright",legend=c("Spp 1", "Spp 2"),lty="solid",col=c("steelblue","firebrick"),bty="n")
}