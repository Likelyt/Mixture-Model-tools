#'pdfmix of the mixture of Cauchy
#'@param x Cauchy mixture data
#'@param mu loaction parameters
#'@param sigma scale parameters
#'@param w weight parameters
#'@return pdfmix of Cauchy mixture data
#'@export
pdfmix=function(x,mu,sigma,w){
  sum(w*dcauchy(x,mu,sigma))
}


#'density of the mixture of Cauchy
#'@param x Cauchy mixture data
#'@param loc loaction parameters
#'@param scale scale parameters
#'@param weigth weight parameters
#'@return density of Cauchy mixture data
#'@export
dmixcauchy=function(x,loc,scale,weight){
  sapply(x,pdfmix,mu=loc,sigma=scale,w=weight)
}


#'cdfmix of the mixture of cauchy
#'@param x Cauchy mixture data
#'@param mu loaction parameters
#'@param sigma scale parameters
#'@param w weight parameters
#'@return cdfmix of Cauchy mixture data
#'@export
cdfmix=function(x,mu,sigma,w){
  sum(w*pcauchy(x,mu,sigma))
}


#'distribution of the mixture of cauchy
#'@param x Cauchy mixture data
#'@param loc loaction parameters
#'@param scale scale parameters
#'@param weight weight parameters
#'@return distribution of Cauchy mixture data
#'@export
pmixcauchy=function(x,loc,scale,weight){
  sapply(x,cdfmix,mu=loc,sigma=scale,w=weight)
}

