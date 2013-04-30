#Stochastic Lotka-Volterra logistic growth model
library(mvtnorm)
library(snow)

cl <- makeCluster(cores, type = "SOCK")

print(date())
# parameters -----------------------------------------
  # name of output file:
  outfile =  "LV_2spp_stoch.csv"

  #specify combinations of freq-dep parameters
  r1_v=c(1,1.5,2)
  r2_v=c(1,1.5,2)

  K1_v=c(50,80,100)
  K2_v=c(50,80,100)

  a12_v=seq(0.3,1,by=0.1)
  a21_v=seq(0.3,1,by=0.1)

  simul0 = expand.grid("r1"=r1_v,"r2"=r2_v,"K1"=K1_v,"K2"=K2_v,"a12"=a12_v,"a21"=a21_v)
  simul = simul0[which(simul[,1]<=simul[,2]),]
  maxtime=5000
  initialNsp1 = 5
  initialNsp2 = 5
  a11=1
  a22=1

  #set up iterations
    iterations=2000
    persistence = matrix(NA,nrow=iterations, ncol=2)   #to save persistence data
    mt = matrix(NA,nrow=iterations, ncol=2)            #to save population data
    saveit = array(NA,dim=c(iterations,2,dim(simul)[1]))
    
#2-species lotka volterra model----------------------------------------------
updateN = function(r_arg, K_arg, Nself, N1, a_intra, a1){
    newN = Nself + Nself * (r_arg * ((K_arg - a_intra*
    Nself - a1 * N1) / K_arg))
    newN=rpois(1,newN)         #demographic stochasticity
    if(newN < 0) newN = 0
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

    r1=simul$r1[iSim]
    r2=simul$r2[iSim]
    K1=simul$K1[iSim]           #get parameter combo
    K2=simul$K2[iSim]
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
    	N[counter,1] = updateN(r[1],K[1],N[i-1,1],N[i-1,2]a11,a12)
    	N[counter,2] = updateN(r[2],K[2],N[i-1,2],N[i-1,1]a22,a21)
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

simul=parsapply(cl, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
write.table(simul,outfile, sep=",",row.names=FALSE)

stopCluster(cl)
