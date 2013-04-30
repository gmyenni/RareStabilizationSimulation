#This script uses saved combinations of frequency-dependent parameters for which
#the equivalence, stabilization, and determinstic results from the Watkinson 
#1980 model have already been calculated (the text file 'simul_rare.txt'). Only 
#parameter combinations in which species 1 was less than or equal to 25% of the 
#total community size in the deterministic solution are used. (Parameter 
#combinations in which lambda1 or lambda2 was less than 15 are run and results 
#are saved, but are not included in the manuscript, due to the confounding 
#effects of low lambdas on extinction times.) 

#This script uses the package snow to run simulations in parallel. The 
#simulations loop is written as a function to accomodate the snow package 
#('write simulation function' below). Within the simulation function, each 
#simulation is iterated 2000 times ('iterations loop' below). Within each 
#iteration, community growth for every parameter combination is simulated until 
#one species goes extinct ('population growth' below, a maximum of 50000 
#timesteps was set as a safeguard, but was never reached). The population growth
#function is written in the file 'updateN_function.r.' Coexistence times and
#mean population sizes are saved for every iteration, and summary statistics are
#calculated at the end ('iteration stats' below).  

print(date())
# parameters -----------------------------------------
  # name of output file:
  outfile =  "annplant_2spp_stoch1.csv"

  #load saved combinations of frequency-dependent parameters
  #in which species one is rare in the deterministic solution
  simul0=read.csv("simul_rare.txt")
  batch=1:1408
  simul=simul0[batch,]         #run current batch
  iterations=2000

#set up simulation output
simul$meanN1 = NA
simul$medN1 = NA
simul$meanN2 = NA
simul$medN2 = NA
simul$sdN1 = NA
simul$sdN2 = NA
simul$MTC = NA
simul$medTC = NA
simul$MTP1 = NA
simul$medTP1 = NA     
simul$MTP2 = NA
simul$medTP2 = NA
simul$Ext1 = NA
simul$Ext2 = NA

#set up parallelization and load population growth function onto 
#each core (see updateN_function.r)
library(snow)
cores=8
cluster1 <- makeCluster(cores, type = "SOCK")
clusterCall(cluster1, function() source("updateN_function.r"))    
                                        

#write simulation function-------------------------------------
Sim <- function(k, simul, iterations) {

    l1=simul$l1[k]           #get parameter combo
    l2=simul$l2[k]
    a11=simul$a11[k]
    a12=simul$a12[k]
    a21=simul$a21[k]
    a22=simul$a22[k]

    #set up iterations
    maxtime=50000
    initialNsp1 = 5
    initialNsp2 = 5
    
    persistence = matrix(NA,nrow=iterations, ncol=2)   #to save persistence data
    mt = matrix(NA,nrow=iterations, ncol=2)            #to save population data

	#iterations loop------------------------------------
	for(j in 1:iterations) {

  	N=data.frame("spp1"=initialNsp1,"spp2"=initialNsp2)  #initialize

		# population growth (see updateN function) ---------------------------
  
   		counter=1
    		stopRun=F
    		while(stopRun==F){
      		counter=counter+1
    			N[counter,1] = updateN(l1,N[counter-1,1],N[counter-1,2],a11,a12) 
    			N[counter,2] = updateN(l2,N[counter-1,2],N[counter-1,1],a22,a21)
      		if(sum(N[counter,]==0)>=1) stopRun=T      #continue until one species goes extinct
      		if(counter>maxtime) stopRun=T             ##or shut down after 50000 timesteps
      			} #next timestep (counter)

  	tmp = rep(dim(N)[1],2)
    	tmp[which(N[dim(N)[1],]>0)]=NA
    	persistence[j,] = tmp     #save persistence data
  	mt[j,] = colMeans(N)        #save mean population sizes
		} #next iteration (j)

#iterations stats
coexist=apply(persistence,MARGIN=1,FUN=min,na.rm=T)
simul$meanN1[k] = mean(mt[,1])                  #mean/median population sizes
simul$medN1[k] = median(mt[,1])
simul$meanN2[k] = mean(mt[,2])
simul$medN2[k] = median(mt[,2])
simul$sdN1[k] = sd(mt[,1])                      #population standard deviations
simul$sdN2[k] = sd(mt[,2])
simul$MTC[k] = mean(coexist,na.rm=T)             #mean/median time populations coexist
simul$medTC[k] = median(coexist,na.rm=T)
simul$MTP1[k] = mean(persistence[,1],na.rm=T)        #mean/median times populations persist
simul$medTP1[k] = median(persistence[,1],na.rm=T)
simul$MTP2[k] = mean(persistence[,2],na.rm=T)
simul$medTP2[k] = median(persistence[,2],na.rm=T)
simul$Ext1[k] = sum(is.na(persistence[,1])==F)/iterations    #fraction of times species goes extinct
simul$Ext2[k] = sum(is.na(persistence[,2])==F)/iterations

#write line
simul[k,]
} #end Sim function

#run simulations in parallel
k=seq(1:dim(simul)[1])
simul1=parSapply(cluster1, k, Sim, simul=simul, iterations=iterations)
simul1=t(simul1)
write.table(simul1,outfile, sep=",",row.names=FALSE)

stopCluster(cluster1)