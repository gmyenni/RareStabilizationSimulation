#2-species annual plant model----------------------------------------------
updateN = function(r_arg, Nself, N1, a_intra, a1){
    newN = (r_arg * Nself) /(1 + a_intra * Nself + a1 * N1 )
    newN=rpois(1,newN)         #demographic stochasticity
    return(newN)
}

  maxtime=50000
  initialNsp1 = 5
  initialNsp2 = 5
  a11=1
  a22=1
  
  l1=12
  l2=14
  a12=0.2
  a21=0.3
  
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

      matplot(N,type="l")