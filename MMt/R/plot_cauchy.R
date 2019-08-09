#'Histogram of Cauchy mixture with true paramter
#'@param x mixture data
#'@param mu lcoation parameters
#'@param sigma scale parameters
#'@param lambda weight parameters
#'@param a plot start from
#'@param b plot end with
#'@return histogram
#'@export
plot_cauchy = function(x, mu, sigma, lambda, a, b){
  temp <- seq(a, b, length.out=1000)
  hist(x, col="lightblue",freq=FALSE, breaks=100, xlim = c(a,b))
  lines(density(x), col="blue", lwd=2) # the empirical density

  l = 0
  for(i in 1:length(lambda)){
    l = l + lambda[i] * dcauchy(temp, mu[i], sigma[i])
  }
  lines(temp, l, col = "black", lty = 2, lwd = 2) #estimate density

  legend("topleft", legend=c("empirical density", "estimate density"),
         col=c("blue", "black"), lty=1:2, cex=0.8)
}
