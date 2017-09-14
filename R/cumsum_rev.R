######function for reverse cusum######
cumsum_rev=function(x){
  y = rev(cumsum(rev(x)))
  return(y)
}
#####################################