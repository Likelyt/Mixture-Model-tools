#'Cauchy mixture data generation
#'@param y mixture data
#'@param m initial guess of mixture components
#'@param general_function general pdf function of components comes from
#'@param threshold relative change ratio of loss
#'@param weight_omit threshold for too small weight, default is 0.01
#'@return number of components, location paramters, scale parameters, and weight paramters
#'@export
NI_QCP = function(y, m, general_function, threshold = 0.001,  weight_omit=0.01){
  require(changepoint)
  require(Rsolnp)
  library(Rsolnp, quietly = T)
  is.integer0 <- function(x){is.integer(x) && length(x) == 0L}

  n = length(y)

  #intialize mu
  mu_hat = rep(0, m)
  mu_hat = as.numeric(quantile(y, seq(1,m,1)/(m+1)))

  #initialize sigma
  sigma_quantile = as.numeric(quantile(y, seq(1,m+1,1)/(m+2)))
  sigma_hat =  rep(0, m)
  for(k in 1:m){
    sigma_hat[k] = (sigma_quantile[k+1] - sigma_quantile[k])/6
  }

  #initialize lambda
  lambda_hat = rep(0,m)
  Fn = ecdf(y) #ecdf
  b = Fn(y)
  A = matrix(0, n, m)
  for(i in 1:n){
    for(j in 1:m){
      A[i,j] = pcauchy(y[i], mu_hat[j], sigma_hat[j])
    }
  }
  lambda_hat = as.vector(lsfit(A, b)$coef[1:m+1])/sum(as.vector(lsfit(A, b)$coef[1:m+1]))

  ansmean=cpt.mean(diff(mu_hat))
  cpt_parameter = min(as.vector(unlist(ansmean@param.est)))
  #plot(ansmean, type = "l", cpt.col = "blue", xlab = "Index",cpt.width = 4)
  cutoff_index = 1

  s = c()
  if(diff(mu_hat)[1] > cpt_parameter){ #begin from positive
    k = 1
    for(i in 1:(m-1)){
      if(diff(mu_hat)[i]<cpt_parameter){
        if(k > 0){
          s[k] = i-1
          k = k * -1
        }
      }
      if(diff(mu_hat)[i]>cpt_parameter){
        if(k < 0){
          s[abs(k)+1] = i-1
          k = k * -1
        }
      }
      i = i + 1
    }
  }else if(diff(mu_hat)[1] < cpt_parameter){
    s = which(diff(mu_hat)>cpt_parameter)-1
  }


  group = length(s)+1
  l_group = vector("list", group)

  for(i in 1:group){
    if(i>1 & i<group){
      l_group[[i]] = (s[i-1]+2):(s[i]+1)
    }else if(i == 1){
      l_group[[i]] = 1:(s[i]+1)
    }else if(i == group){
      l_group[[i]] = (s[i-1]+2):m
    }
  }

  mean_mu =rep(0,group)
  mean_sigma = rep(0,group)
  mean_lambda = rep(0,group)
  for(i in 1:group){
    mean_mu[i] = mean(mu_hat[l_group[[i]]])
    mean_sigma[i] = mean(sigma_hat[l_group[[i]]])
    #mean_sigma[i] = sqrt(sum(sigma_hat[l_group[[i]]]^2))
    mean_lambda[i] = sum(lambda_hat[l_group[[i]]])
  }

  #scale weight
  mean_lambda = mean_lambda/sum(mean_lambda)

  #step 2 Updtae parameters
  K = group

  log_li_mix_sigma = function(sigma, mu, weight){
    l = 0
    for(i in 1:n){l = l + log(sum(weight*dcauchy(y[i], mu, sigma)))}
    return(-l)
  }

  log_li_mix_mu = function(mu, sigma, weight){
    l = 0
    for(i in 1:n){l = l + log(sum(weight*dcauchy(y[i], mu, sigma)))}
    return(-l)
  }

  log_li_mix_weight = function(weight, mu, sigma){
    l = 0
    for(i in 1:n){l = l + log(sum(weight*dcauchy(y[i], mu, sigma)))}
    return(-l)
  }

  fheq = function(x, mu, sigma) {sum(x)}


  log_mix = function(mu, sigma, weight){
    l = 0
    for(i in 1:n){l = l + log(sum(weight*dcauchy(y[i], mu, sigma)))}
    return(-l)
  }


  #step 2.2 Initial values
  New_mu = matrix(0,nrow = 100, ncol = K)
  New_sigma = matrix(0,nrow = 100, ncol = K)
  New_weight = matrix(0,nrow = 100, ncol = K)

  New_mu[1,] = mean_mu
  New_sigma[1,] = mean_sigma
  New_weight[1,] = mean_lambda

  STOP = 100
  for(i in 2:100){
    #For i update sigma
    New_sigma[i,] = optim(par = rep(1,K), fn = log_li_mix_sigma, method = "L-BFGS-B",
                          lower = rep(0.0001,K), mu = New_mu[i-1,],
                          weight = New_weight[i-1,])$par

    #For i update mu
    New_mu[i,] = optim(par = New_mu[i-1,], fn = log_li_mix_mu, method = "L-BFGS-B",
                       sigma = New_sigma[i,], weight = New_weight[i-1,])$par

    #For i update weight
    New_weight[i,] = solnp(New_weight[i-1,], log_li_mix_weight,
                           eqfun = fheq, eqB = 1,
                           LB = rep(0,K),UB = rep(1,K),
                           mu=New_mu[i,], sigma = New_sigma[i,], control = list(trace = 0))$pars

    #Check the difference
    if(abs(log_mix(New_mu[i,],New_sigma[i,], New_weight[i,])
           - log_mix(New_mu[i-1,],New_sigma[i-1,], New_weight[i-1,])) <
       threshold * log_mix(New_mu[i-1,],New_sigma[i-1,], New_weight[i-1,])){
      STOP = i
      break
    }
  }

  if(!is.integer0(which(New_weight[STOP,] < weight_omit))){
    omit_comp = which(New_weight[STOP,] < weight_omit)
    group = length(New_weight[STOP,]) - length(omit_comp)
    #rescale weight
    final_sigma = New_sigma[STOP,][-omit_comp]
    final_mu = New_mu[STOP,][-omit_comp]
    final_weight = New_weight[STOP,][-omit_comp]
    final_weight = final_weight/sum(final_weight)
    return(list(k = group, mu = final_mu, sigma = final_sigma, lambda = final_weight))
  }else{
    return(list(k = group, mu = New_mu[STOP,], sigma = New_sigma[STOP,], lambda = New_weight[STOP,]))
  }

}



