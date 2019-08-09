#'Degree of overlapping (DOL) function for mixture model
#'@param x Cauchy mixture data
#'@param alpha loaction parameters
#'@param beta scale parameters
#'@param w weight parameters
#'@return DOL function
#'@export
DOL_f <- function(x, alpha, beta, w){
  k = length(alpha)
  f = rep(NA, k)
  for(i in 1:k){
    assign(paste0('f',i,sep=""), w[i]*dcauchy(x, alpha[i], beta[i]))
  }
  temp = 10^7
  for(i in 1:k){
    temp = pmin(temp, eval(parse(text=paste0('f',i,sep=""))))
  }
  return(temp)
}

#'Degree of overlapping (DOL)  for mixture model
#'@param alpha loaction parameters
#'@param beta scale parameters
#'@param w weight parameters
#'@return DOL
#'@export
DOL <-function(alpha, beta, w){
  dol = as.numeric(integrate(DOL_f, -Inf, Inf, alpha, beta, w)[1])/min(w)
  return(dol)
}


#'A measure of between components dispersion (BCD) for mixture model
#'@param alpha loaction parameters
#'@param beta scale parameters
#'@param w weight parameters
#'@return BCD
#'@export
BCD <- function(alpha, beta, w){
  set.seed(12345)
  n = 10^5
  k = length(alpha)
  y = rmixcauchy(n, alpha, beta, w, seed = 12345)
  num = 0
  for(i in 1:k){
    num = num + w[i]*abs(alpha[i]-median(y))
  }
  den = (1/n)*sum(abs(y - median(y)))
  bcd = num/den
  return(bcd)
}


