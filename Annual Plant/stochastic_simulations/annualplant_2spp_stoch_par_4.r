library(snow)
cores=8

chunk=1:792
#chunk=793:1584
#chunk=1585:2376

print(date())
# parameters -----------------------------------------
  # name of output file:
  outfile =  "annplant_2spp_stoch_4.csv"

  #specify combinations of freq-dep parameters
  l1_v=10:20
  l2_v=10:20

  a12_v=seq(0.1,1,by=0.1)
  a21_v=seq(0.1,1,by=0.1)

  simul0 = expand.grid("l1"=l1_v,"l2"=l2_v,"a12"=a12_v,"a21"=a21_v)
  simul = simul0[which(simul0[,1]<=simul0[,2]),]

  lista21=which(simul$a21<=0.2)
  lista12=which(simul$a12<=0.2)
  listall=unique(c(lista21,lista12))
  simul = simul[listall,]
  simul = simul[chunk,]
  
  iterations=2000
    
#simulation output
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

cluster1 <- makeCluster(cores, type = "SOCK")
clusterCall(cluster1, function() source("updateN_function.r"))

#simulation function-------------------------------------
Sim <- function(k, simul, iterations) {

    l1=simul$l1[k]           #get parameter combo
    l2=simul$l2[k]
    a12=simul$a12[k]
    a21=simul$a21[k]
    a11=1
    a22=1

    #set up iterations
    maxtime=10000
    initialNsp1 = 5
    initialNsp2 = 5
    
    persistence = matrix(NA,nrow=iterations, ncol=2)   #to save persistence data
    mt = matrix(NA,nrow=iterations, ncol=2)            #to save population data

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
      if(counter>maxtime) stopRun=T
      } #next timestep (counter)

  tmp = rep(dim(N)[1],2)
    tmp[which(N[dim(N)[1],]>0)]=NA
    persistence[j,] = tmp     #save persistence data
  mt[j,] = colMeans(N)        #save mean population sizes
} #next iteration (j)

#iterations stats
#simul=c(l1,l2,a12,a21)
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

#progress update
  #tmp = paste("Parameter combination",k,"of",dim(simul)[1],sep=" ")
  #print(tmp)
  #print(date())
  #flush.console()
  
#write line
  simul[k,]
} #end Sim function

k=seq(1:dim(simul)[1])
simul1=parSapply(cluster1, k, Sim, simul=simul, iterations=iterations)
simul1=t(simul1)
write.table(simul1,outfile, sep=",",row.names=FALSE)

stopCluster(cluster1)