library(snow)
cores=8
print(date())
# parameters -----------------------------------------
  # name of output file:
  outfile =  "annplant_2spp_det1.csv"

  #specify combinations of freq-dep parameters
  l1_v=10:20
  l2_v=10:20

  a11_v=c(0.1,0.3,0.5)
  a12_v=c(0.1,0.3,0.5,0.7,0.9,1)
  a21_v=c(0.1,0.3,0.5,0.7,0.9,1)
  a22_v=c(0.1,0.3,0.5,0.7,0.9,1)

  simul = expand.grid("l1"=l1_v,"l2"=l2_v,"a11"=a11_v,"a12"=a12_v,"a21"=a21_v,"a22"=a22_v)

#simulation output
simul$N1 = NA
simul$N2 = NA

cluster1 <- makeCluster(cores, type = "SOCK")
clusterCall(cluster1, function() source("analyN_function.r"))

#simulation function-------------------------------------
Sim = function(k, simul) {

    l1=simul$l1[k]           #get parameter combo
    l2=simul$l2[k]
    a11=simul$a11[k]
    a12=simul$a12[k]
    a21=simul$a21[k]
    a22=simul$a22[k]
    
#simulation stats
simul$N1[k] = analyN(l1,l2,a11,a12,a21,a22)[1]
simul$N2[k] = analyN(l1,l2,a11,a12,a21,a22)[2]

simul[k,]
}

k=seq(1:dim(simul)[1])
simul1=parSapply(cluster1, k, Sim, simul=simul)
simul1=t(simul1)
write.table(simul1,outfile, sep=",",row.names=FALSE)
stopCluster(cluster1)
