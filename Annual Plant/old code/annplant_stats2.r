 #specify combinations of freq-dep parameters

#frequency dependence
l1_v=16:18
l2_v=13:15
l3_v=10:12

a12_v=seq(0.1,0.3,by=.05)
a13_v=seq(0.1,0.3,by=.05)
a21_v=seq(0.3,0.4,by=.05)
a23_v=seq(0.3,0.4,by=.05)
a31_v=seq(0.4,0.6,by=.05)
a32_v=seq(0.4,0.6,by=.05)

annplant_stat2=matrix(NA, nrow=1, ncol=15)

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

    source("annualplant2.r")

annplant_stat2=rbind(annplant_stat2,c(l1,l2,l3,a12,a13,a21,a23,a31,a32,mean(N[,1]),mean(N[,2]),mean(N[,1]),max(which(N[,1] > 0)),max(which(N[,2] > 0)),max(which(N[,3] > 0))))
write.csv(annplant_stat2,file="annplant_stat3")
}}}}}}}}}


