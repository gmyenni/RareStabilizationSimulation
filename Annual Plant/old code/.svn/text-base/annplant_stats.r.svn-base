#specify combinations of freq-dep parameters

#frequency dependence
l1_v=15:20
l2_v=13:17
l3_v=10:15

a12_v=seq(0.1,0.9,by=.2)
a13_v=seq(0.1,0.9,by=.2)
a21_v=seq(0.1,0.9,by=.2)
a23_v=seq(0.1,0.9,by=.2)
a31_v=seq(0.1,0.9,by=.2)
a32_v=seq(0.1,0.9,by=.2)

annplant_stat=matrix(NA, nrow=1, ncol=24)

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

    source("annplant_iterations.r")

annplant_stat=rbind(annplant_stat,c(l1,l2,l3,a12,a13,a21,a23,a31,a32,mean(mt[,1]),mean(mt[,2]),mean(mt[,3]),sd(mt[,1]),sd(mt[,2]),sd(mt[,3]),pop1MTE,pop2MTE,pop3MTE,wins1,wins2,wins3))

write.csv(annplant_stat,file="annplant_stat")

}}}}}}}}}

