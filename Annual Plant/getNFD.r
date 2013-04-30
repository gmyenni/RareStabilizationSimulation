
# This is a function to calculate negative frequency dependence
# for each species in a two species annual plant  model.
# You input the parameters (r's, alpha's) and it returns
# the growth rate of each species at very low density,
# the frequency of each species at low density,
# the growth rate at a slightly higher density,
# and the frequency at the higher density.

# Main note: When I increase the growth rate of species 1 
# a little above low density, I set species 2 at its equilibrium
# density, given the "fixed" density of species 1.

# For comparison, I also included the function we used in the paper.


# FUNCTIONS ----------------------------------------------------

# strength of stabilization as measured in our paper
SOS=function(r,alpha){ 
  out=rep(NA,2)
  out[1]=r[2]/(1+(alpha[1,2]/alpha[2,2])*(r[2]-1))
  out[2]=r[1]/(1+(alpha[2,1]/alpha[1,1])*(r[1]-1))
  return(out)
}


# calculate log per capita growth rate
getPCG=function(r,alpha,N){
  # r is a vectory[2] of fecundities for spp 1 and 2
  # alpha is a matrix[2,2] of interaction coefficients where alpha[i,j] is the effect of j on i
  # N is a vectory[2] of densities for spp 1 and 2
  newN=rep(NA,2)
  newN[1] <- r[1]*N[1] / (1 + alpha[1,1]*N[1]+alpha[1,2]*N[2])
  newN[2]<- r[2]*N[2] / (1 + alpha[2,1]*N[1]+alpha[2,2]*N[2])
  out=log(newN)-log(N) # vector[2] of per capita growth rates
  return(out)
}

# figure out the equilibrium density of the focal species' competitor
# given a fixed density of the focal species
getEqDensity=function(species,r,alpha,N.star){
  # species = focal species (this one is fixed, density of the other is returned)
  # N.start (a scalar) is the density of the fixed spp
	if(species==1)  
	{	# density of species one is known
		out<-(r[2]-1-alpha[2,1]*N.star)/alpha[2,2]
	}else{
		# density of species two is known
		out<-(r[1]-1-alpha[1,2]*N.star)/alpha[1,1]
	}
	return(out)  
}

getNFD=function(r,alpha,lowN,deltaN){
 # low N is the low density 
 # and deltaN is the increase in N for the higher density to explore
 
 # vectors for output
  pgr1=freq1=pgr2=freq2=rep(NA,2)
 
 # get low density growth rate for spp 1
 tmpN=rep(NA,2)
 tmpN[1]=lowN  # set density of focal spp
 tmpN[2]=getEqDensity(species=1,r,alpha,N.star=tmpN[1])  # figure out density of its competitor
 tmpOut=getPCG(r,alpha,N=tmpN)
 pgr1[1]=tmpOut[1]
 freq1[1]=tmpN[1]/sum(tmpN)
  
 # get low density growth rate for spp 2
 tmpN=rep(NA,2)
 tmpN[2]=lowN
 tmpN[1]=getEqDensity(species=2,r,alpha,N.star=tmpN[2])
 tmpOut=getPCG(r,alpha,N=tmpN)
 pgr1[2]=tmpOut[2]
 freq1[2]=tmpN[2]/sum(tmpN) 
  
# get higher density growth rate for spp1 1
 tmpN=rep(NA,2)
 tmpN[1]=lowN + deltaN
 tmpN[2]=getEqDensity(species=1,r,alpha,N.star=tmpN[1])
 tmpOut=getPCG(r,alpha,N=tmpN)
 pgr2[1]=tmpOut[1]
 freq2[1]=tmpN[1]/sum(tmpN)
  
 # get higher density growth rate for spp1 2
 tmpN=rep(NA,2)
 tmpN[2]=lowN + deltaN
 tmpN[1]=getEqDensity(species=2,r,alpha,N.star=tmpN[2])
 tmpOut=getPCG(r,alpha,N=tmpN)
 pgr2[2]=tmpOut[2]
 freq2[2]=tmpN[2]/sum(tmpN)
  
 return(list(pgr1=pgr1,freq1=freq1,pgr2=pgr2,freq2=freq2))
}


# example use --------------------------------------------------------

r=c(100,80)
alpha=cbind(c(1,0.2),c(0.8,1))
lowN=0.001
deltaN=10  # results not very sensitive to deltaN until it gets large--suggests some nonlinearity in NFD

# get the frequencies and growth rates
test=getNFD(r,alpha,lowN,deltaN)

# calculate the NFD
# NFD = rise over run * -1 
print(-1*(test$pgr2-test$pgr1)/(test$freq2-test$freq1))  

# calculate Yenni et al. 2012 SOS
print(SOS(r,alpha)) # compare to what we used in the paper

# take a look
matplot(x=rbind(test$freq1,test$freq2),y=rbind(test$pgr1,test$pgr2),type="b",
      xlab="Frequency",ylab="PGR")

# compare SOS and NFD for 1000 random parameter sets
tmp=rep(NA,1000)
compareRandom=data.frame(SOS1=tmp,SOS2=tmp,NFD1=tmp,NFD2=tmp)
for(i in 1:1000){
  r=rnorm(2,100,10)
  alpha=rnorm(4,1,0.2)
  alpha=matrix(alpha,2,2)
  compareRandom[i,1:2]=SOS(r,alpha)
  test=getNFD(r,alpha,lowN,deltaN)
  compareRandom[i,3:4]=-1*(test$pgr2-test$pgr1)/(test$freq2-test$freq1)
}
par(mfrow=c(1,2))
plot(compareRandom$SOS1,compareRandom$NFD1,xlab="SOS",ylab="NFD",main="Spp 1")
plot(compareRandom$SOS2,compareRandom$NFD2,xlab="SOS",ylab="NFD",main="Spp 2")
