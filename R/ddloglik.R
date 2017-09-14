#######function for calculate L1,L2#######
ddloglik= function(n,delta,z,beta){
  p = length(beta)
  
  scaled_schoenfeld=matrix(0,nrow=n,ncol=p)
  #GVG=matrix(0,nrow=p*knot,ncol=p*knot)
  S0=rev(cumsum(rev(exp(z%*%beta))))
  S1_pre=sweep(z, MARGIN=1, exp(z%*%beta),'*')
  S1=apply(S1_pre,MARGIN=2,cumsum_rev)
  S2=rep(0,n*p*p)
  dim(S2)=c(n,p,p)
  #B_square=rep(0,n*knot*knot)
  #dim(B_square)=c(n,knot,knot)
  V=S2
  L1=colSums(delta*(z - S1/S0))
  partial_likelihood=delta*((z%*%beta)-log(S0))
  
  L2=matrix(0, nrow=p,ncol=p)
  j=0
  repeat{
    j=j+1
    S2[,,j]=apply(matrix(z*as.vector(exp(z%*%beta))*z[,j],ncol=p),MARGIN=2,cumsum_rev)
    V[,,j]=S2[,,j]/S0-S1*S1[,j]/(S0**2)
    L2[,j]=matrix(colSums(delta*V[,,j]),ncol=1)
    if (j==p) break
  }
  Lambda=cumsum(delta/S0)
  list(L1=L1,L2=L2, partial_likelihood=partial_likelihood,Lambda=Lambda,S2=S2) 
} 