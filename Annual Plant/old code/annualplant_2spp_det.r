print(date())
# parameters -----------------------------------------
  # name of output file:
  outfile =  "annplant_2spp_det.csv"

  #specify combinations of freq-dep parameters
  l1_v=10:20
  l2_v=10:20

  a12_v=seq(0,1,by=0.1)
  a21_v=seq(0,1,by=0.1)

  simul = expand.grid("l1"=l1_v,"l2"=l2_v,"a12"=a12_v,"a21"=a21_v)
  maxtime=500
  initialNsp1 = 5
  initialNsp2 = 5
  a11=1
  a22=1

#2-species annual plant model----------------------------------------------
updateN = function(r_arg, Nself, N1, a_intra, a1){
    newN = (r_arg * Nself) /(1 + a_intra * Nself + a1 * N1 )
    return(newN)
}

#simulation output
simul$meanN1 = NA
simul$medN1 = NA
simul$meanN2 = NA
simul$medN2 = NA
simul$sdN1 = NA
simul$sdN2 = NA
simul$coexist = NA
simul$persist1 = NA
simul$persist2 = NA

#simulation loop-------------------------------------
for(iSim in 1:dim(simul)[1]) {

    l1=simul$l1[iSim]           #get parameter combo
    l2=simul$l2[iSim]
    a12=simul$a12[iSim]
    a21=simul$a21[iSim]

  N=data.frame("spp1"=initialNsp1,"spp2"=initialNsp2)  #initialize

# populations loop ---------------------------------
    counter=1
    stopRun=F
    while(stopRun==F){
      counter=counter+1
    	N[counter,1] = updateN(l1,N[counter-1,1],N[counter-1,2],a11,a12)
    	N[counter,2] = updateN(l2,N[counter-1,2],N[counter-1,1],a22,a21)
      if(sum(N[counter,]<1)>=1) stopRun=T
      if(counter==maxtime) stopRun=T
      } #next timestep (counter)

#simulation stats
simul$meanN1[iSim] = mean(N[,1])
simul$medN1[iSim] = median(N[,1])
simul$meanN2[iSim] = mean(N[,2])
simul$medN2[iSim] = median(N[,2])
simul$sdN1[iSim] = sd(N[,1])
simul$sdN2[iSim] = sd(N[,2])
simul$coexist[iSim] = counter
simul$persist1[iSim] = floor(N[counter,1])
simul$persist2[iSim] = floor(N[counter,2])

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