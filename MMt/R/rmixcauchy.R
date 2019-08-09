#'Cauchy mixture data generation
#'@param n number of sample
#'@param loc loaction parameters
#'@param scale scale parameters
#'@param weigth weight parameters
#'@param seed random seed
#'@return Cauchy mixture data
#'@export
rmixcauchy = function(n, loc, scale, weight, seed){
  set.seed(seed)
  m = length(weight)
  z = sample.int(m,n,replace=T,prob=weight)
  x= rcauchy(n,loc[z],scale[z])
  return(x)
}
