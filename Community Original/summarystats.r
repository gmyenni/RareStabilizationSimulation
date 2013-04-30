#specify combinations of freq-dep parameters

#frequency dependence
sp1_v=-1:-3                  #species 1 slope
sp2_v=-1:-5
sp3_v=-1:-7
int1_v=seq(.4,.6,by=.1)      #species 1 intercept
int2_v=seq(.3,.5,by=.1)
int3_v=seq(.1,.4,by=.1)

stat=matrix(NA, nrow=1, ncol=15)

for(a in 1:length(sp1_v)) {
  for(b in 1:length(int1_v)) {
  for(c in 1:length(sp2_v)) {
  for(d in 1:length(int2_v)) {
  for(e in 1:length(sp3_v)) {
  for(f in 1:length(int3_v)) {
  
    sp1=sp1_v[a]
    int1=int1_v[b]
    sp2=sp2_v[c]
    int2=int2_v[d]
    sp3=sp3_v[e]
    int3=int3_v[f]
    
    source("iterations.r")
  
stat=rbind(stat,c(sp1,int1,sp2,int2,sp3,int3,mean(mt[,1]),mean(mt[,2]),mean(mt[,3]),sd(mt[,1]),sd(mt[,2]),sd(mt[,3]),wins1,wins2,wins3))

}}}}}}

write.csv(stat,file="stat")