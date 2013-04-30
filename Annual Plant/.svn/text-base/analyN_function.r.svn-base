#2-species annual plant model analytical solution----------------------------------------------
analyN = function(r1, r2, a1, a12, a21, a2){
    N1 = (r1-1-(a12/a2)*(r2-1))/(a1-a21*a12/a2)
    N2 = (r2-1-(a21/a1)*(r1-1))/(a2-a21*a12/a1)


if(sum(is.infinite(N1),is.infinite(N2),is.nan(N1),is.nan(N2))>0) {

initialNsp1=0
initialNsp2=0

N=data.frame("N1"=initialNsp1,"N2"=initialNsp2)  #initialize

for(i in 2:100) {

N[i,1] = max((r1-1-a12*N[i-1,2])/a1,0)
N[i,2] = max((r2-1-a21*N[i-1,1])/a2,0)

N1=apply(N,2,mean)[1]
N2=apply(N,2,mean)[2]

  } }
  
if(N1<0) { N1 = 0; N2 = (r2-1)/a2  }

if(N2<0) { N2 = 0; N1 = (r1-1)/a1  }
    
    return(c(N1,N2))
}