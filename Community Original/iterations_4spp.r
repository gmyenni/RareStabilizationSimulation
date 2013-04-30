#NOTE: this code requires that the two scripts be saved in the same directory.

iterations=200
xt=matrix(NA, nrow=iterations, ncol=4)
mt=matrix(NA,nrow=iterations, ncol=4)


for(j in 1:iterations) {
  source("stochastic_freq-depend_4spp.r")
  xt[j,1] = max(which(N[,1] > 0))
  xt[j,2] = max(which(N[,2] > 0))
  xt[j,3] = max(which(N[,3] > 0))
  xt[j,4] = max(which(N[,4] > 0))
  mt[j,1] = mean(N[,1])
  mt[j,2] = mean(N[,2])
  mt[j,3] = mean(N[,3])
  mt[j,4] = mean(N[,4])
}

#mean and standard deviation of time to extinction for each population
#(filter out all time-to-extinction values that equal the total run-time)
pop1MTE = mean(xt[which(xt[,1] != TotTime),1])
pop2MTE = mean(xt[which(xt[,2] != TotTime),2])
pop3MTE = mean(xt[which(xt[,3] != TotTime),3])
pop4MTE = mean(xt[which(xt[,4] != TotTime),4])

pop1SDTE = sd(xt[which(xt[,1] != TotTime),1])
pop2SDTE = sd(xt[which(xt[,2] != TotTime),2])
pop3SDTE = sd(xt[which(xt[,3] != TotTime),3])
pop4SDTE = sd(xt[which(xt[,4] != TotTime),4])

#number of times each species survives to the end of the simulation
wins1 = length(which(xt[,1] == TotTime))
wins2 = length(which(xt[,2] == TotTime))
wins3 = length(which(xt[,3] == TotTime))
wins4 = length(which(xt[,4] == TotTime))

expected=c(1/4,1/4,1/4,1/4)
wins=c(wins1,wins2,wins3,wins4)

#mean and sd times to extinction for each population.
pop1MTE
pop2MTE
pop3MTE
pop4MTE
pop1SDTE
pop2SDTE
pop3SDTE
pop4SDTE

#print mean and SD population abundance
mean(mt[,1])  #mean pop 1
mean(mt[,2])  #mean pop 2
mean(mt[,3])  #mean pop 3
mean(mt[,4])  #mean pop 4
sd(mt[,1])    #sd pop 1
sd(mt[,2])    #sd pop 2
sd(mt[,3])    #sd pop 3
sd(mt[,4])    #sd pop 4

#no. simulations:
iterations
#number of times that Pop1 survives to the end of the simulation
wins1
#number of times that Pop2 survives to the end of the simulation
wins2
#number of times that Pop2 survives to the end of the simulation
wins3
#number of times that Pop2 survives to the end of the simulation
wins4

#chisq.test(wins,p=expected)