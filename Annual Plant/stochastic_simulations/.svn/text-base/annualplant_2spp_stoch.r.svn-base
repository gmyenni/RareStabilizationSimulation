print(date())
# parameters -----------------------------------------
  # name of output file:
  outfile =  "annplant_2spp_stoch.csv"

  #specify combinations of freq-dep parameters
  l1_v=10
  l2_v=11

  a12_v=0.2
  a21_v=0.3

  simul0 = expand.grid("l1"=l1_v,"l2"=l2_v,"a12"=a12_v,"a21"=a21_v)
  simul = simul0[which(simul0[,1]<=simul0[,2]),]
  maxtime=50000
  initialNsp1 = 5
  initialNsp2 = 5
  a11=1                    
  a22=1

  #set up iterations
    iterations=1
    persistence = matrix(NA,nrow=iterations, ncol=2)   #to save persistence data
    mt = matrix(NA,nrow=iterations, ncol=2)            #to save population data
    saveit = array(NA,dim=c(iterations,2,dim(simul)[1]))
#2-species annual plant model----------------------------------------------
updateN = function(r_arg, Nself, N1, a_intra, a1){
    newN = (r_arg * Nself) /(1 + a_intra * Nself + a1 * N1 )
    newN=rpois(1,newN)         #demographic stochasticity
    return(newN)
}

#simulation output
simul$meanN1 = NA     #mean/median population sizes       
simul$medN1 = NA
simul$meanN2 = NA
simul$medN2 = NA
simul$sdN1 = NA       #population standard deviations        
simul$sdN2 = NA
simul$MTC = NA        #mean/median time populations coexist
simul$medTC = NA
simul$MTP1 = NA       #mean/median times populations persist
simul$medTP1 = NA     
simul$MTP2 = NA
simul$medTP2 = NA
simul$Ext1 = NA       #fraction of times species goes extinct
simul$Ext2 = NA

#simulation loop-------------------------------------
for(iSim in 1:dim(simul)[1]) {

    l1=simul$l1[iSim]           #get parameter combo
    l2=simul$l2[iSim]
    a12=simul$a12[iSim]
    a21=simul$a21[iSim]

#iterations loop------------------------------------
for(j in 1:iterations) {

  N=data.frame("spp1"=initialNsp1,"spp2"=initialNsp2)  #initialize

# populations loop ---------------------------------
    counter=1
    stopRun=F
    while(stopRun==F){
      counter=counter+1
    	N[counter,1] = updateN(l1,N[counter-1,1],N[counter-1,2],a11,a12) 
    	N[counter,2] = updateN(l2,N[counter-1,2],N[counter-1,1],a22,a21)
      if(sum(N[counter,]==0)>=1) stopRun=T
     
      } #next timestep (counter)

  tmp = rep(dim(N)[1],2)
    tmp[which(N[dim(N)[1],]>0)]=NA
    persistence[j,] = tmp     #save persistence data
  mt[j,] = colMeans(N)        #save mean population sizes
} #next iteration (j)

#iterations stats
saveit[,,iSim]=persistence
coexist=apply(persistence,MARGIN=1,FUN=min,na.rm=T)
simul$meanN1[iSim] = mean(mt[,1])             
simul$medN1[iSim] = median(mt[,1])
simul$meanN2[iSim] = mean(mt[,2])
simul$medN2[iSim] = median(mt[,2])
simul$sdN1[iSim] = sd(mt[,1])               
simul$sdN2[iSim] = sd(mt[,2])
simul$MTC[iSim] = mean(coexist,na.rm=T)
simul$medTC[iSim] = median(coexist,na.rm=T)
simul$MTP1[iSim] = mean(persistence[,1],na.rm=T)
simul$medTP1[iSim] = median(persistence[,1],na.rm=T)
simul$MTP2[iSim] = mean(persistence[,2],na.rm=T)
simul$medTP2[iSim]=median(persistence[,2],na.rm=T)
simul$Ext1[iSim]=sum(is.na(persistence[,1])==F)/iterations
simul$Ext2[iSim]=sum(is.na(persistence[,2])==F)/iterations   

#progress update
  tmp = paste("Parameter combination",iSim,"of",dim(simul)[1],sep=" ")
  print(tmp)
  print(date())
  flush.console()
} #next parm combo (iSim)

write.table(simul,outfile, sep=",",row.names=FALSE)


## if only one run performed, make figures
if(dim(simul)[1]==1 & iterations==1){
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