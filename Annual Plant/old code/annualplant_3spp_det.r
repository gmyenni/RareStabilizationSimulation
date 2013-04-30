#specify combinations of freq-dep parameters
l1_v=13:20
l2_v=13:20
l3_v=13:20

a12_v=seq(0.1,0.9,by=.1)
a13_v=seq(0.1,0.9,by=.1)
a21_v=seq(0.1,0.9,by=.1)
a23_v=seq(0.1,0.9,by=.1)
a31_v=seq(0.1,0.9,by=.1)
a32_v=seq(0.1,0.9,by=.1)

#set up
TotTime=200
initialNsp1 = 5
initialNsp2 = 5
initialNsp3 = 5
a11=1
a22=1
a33=1

#3-species annual plant model
updateN = function(r_arg, Nself, N1, N2, a_intra, a1, a2){
    newN = r_arg * Nself /(1 + a_intra * Nself + a1 * N1 + a2 * N2)
    if(newN < 0) newN = 0
    return(newN)
}

#simulation loops-------------------------------------
for(a in 1:length(l1_v)) {
  for(b in 1:length(l2_v)) {
  for(c in 1:length(l3_v)) {
  for(d in 1:length(a12_v)) {
  for(e in 1:length(a13_v)) {
  for(f in 1:length(a21_v)) {
  for(g in 1:length(a23_v)) {
  for(h in 1:length(a31_v)) {
  for(i in 1:length(a32_v)) {

    l1=l1_v[a]
    l2=l2_v[b]
    l3=l3_v[c]

    a12=a12_v[d]
    a13=a13_v[e]
    a21=a21_v[f]
    a23=a23_v[g]
    a31=a31_v[h]
    a32=a32_v[i]

N = matrix(data=NA,nrow=TotTime,ncol=3)
  N[1,1] = initialNsp1
  N[1,2] = initialNsp2
  N[1,3] = initialNsp3

# populations loop ---------------------------------
for(t in 2:TotTime) {
    N[t,1] = floor(updateN(l1,N[t-1,1],N[t-1,2],N[t-1,3],a11,a12,a13))
    N[t,2] = floor(updateN(l2,N[t-1,2],N[t-1,1],N[t-1,3],a22,a21,a23))
    N[t,3] = floor(updateN(l3,N[t-1,3],N[t-1,1],N[t-1,2],a33,a31,a32))
}

write.table(cbind(l1,l2,l3,a12,a13,a21,a23,a31,a32,mean(N[,1]),mean(N[,2]),mean(N[,3]),max(which(N[,1] > 0)),max(which(N[,2] > 0)),max(which(N[,3] > 0))), "annplant_3spp_det.csv", sep=",", append=TRUE, col.names=FALSE)
}}}}}}}}}


