print(date())
library(snow)
cores=8

# parameters -----------------------------------------
  # name of output file:
  outfile =  "annplant_2spp_det.csv"

  #specify combinations of freq-dep parameters
  l1_v=10:20
  l2_v=10:20

  a12_v=seq(0,1,by=0.1)
  a21_v=seq(0,1,by=0.1)

  simul0 = expand.grid("l1"=l1_v,"l2"=l2_v,"a12"=a12_v,"a21"=a21_v)
  simul = simul0[which(simul0[,1]<=simul0[,2]),]
  
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

cluster1 <- makeCluster(cores, type = "SOCK")
clusterCall(cluster1, function() source("updateN_det_function.r"))

#simulation function-------------------------------------
Sim <- function(k, simul, iterations) {

    l1=simul$l1[k]           #get parameter combo
    l2=simul$l2[k]
    a12=simul$a12[k]
    a21=simul$a21[k]
    a11=1
    a22=1
    
  maxtime=50000
  initialNsp1 = 5
  initialNsp2 = 5

  N=data.frame("spp1"=initialNsp1,"spp2"=initialNsp2)  #initialize

# populations loop ---------------------------------
    counter=1
    stopRun=F
    while(stopRun==F){
      counter=counter+1
    	N[counter,1] = updateN(l1,N[counter-1,1],N[counter-1,2],a11,a12)
    	N[counter,2] = updateN(l2,N[counter-1,2],N[counter-1,1],a22,a21)
      if(sum(N[counter,]<1)>=1) stopRun=T
      if(counter>maxtime) stopRun=T
      } #next timestep (counter)

#simulation stats
simul$meanN1[k] = mean(N[,1])
simul$medN1[k] = median(N[,1])
simul$meanN2[k] = mean(N[,2])
simul$medN2[k] = median(N[,2])
simul$sdN1[k] = sd(N[,1])
simul$sdN2[k] = sd(N[,2])
simul$coexist[k] = counter
simul$persist1[k] = floor(N[counter,1])
simul$persist2[k] = floor(N[counter,2])

#progress update
  tmp = paste("Parameter combination",iSim,"of",dim(simul)[1],sep=" ")
  print(tmp)
  print(date())
  flush.console()
  simul[k,]
} #end Sim function


k=seq(1:dim(simul)[1])
simul1=parSapply(cluster1, k, Sim, simul=simul, iterations=iterations)
simul1=t(simul1)
write.table(simul1,outfile, sep=",",row.names=FALSE)

stopCluster(cluster1)
