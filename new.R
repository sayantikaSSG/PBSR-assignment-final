MyMLE = function(data_obtained,modelnumber,initial)
{
  if(modelnumber==1)
  {
    NegLogLike = function(initial,data_obtained)
    {
      alpha = exp(initial[1])
      sigma1 = exp(initial[2])
      mu = initial[3]
      sigma2 = exp(initial[4])
      p = exp(initial[5])/(1 + exp(initial[5]))
      n = length(data_obtained)
      logl = 0
      for(i in 1:n)
      { 
        logl = logl + log(p*dgamma(data_obtained[i],shape = alpha,
                                   scale = sigma1) + (1-p)*dnorm(data_obtained[i],
                                                                 mean = mu,sd=sigma2))
      }
      return(-logl)
    }
    fit = optim(initial,NegLogLike,data_obtained=data_obtained)
    alphar = exp(fit$par[1])
    sigma1r = exp(fit$par[2])
    mur = fit$par[3]
    sigma2r = exp(fit$par[4])
    pr = exp(fit$par[5])/(1 + exp(fit$par[5]))
    toreturn = c(length(fit$par),alphar,sigma1r,mur,sigma2r,pr,-fit$value)
    return(toreturn)
  }
  else if(modelnumber==2)
  {
    NegLogLike = function(initial,data_obtained)
    {
      alpha1 = exp(initial[1])
      sigma1 = exp(initial[2])
      alpha2 = exp(initial[3])
      sigma2 = exp(initial[4])
      p = exp(initial[5])/(1 + exp(initial[5]))
      n = length(data_obtained)
      logl = sum(log(p*dgamma(data_obtained,shape = alpha1,
                              scale = sigma1) + (1-p)*dgamma(data_obtained,
                                                             shape = alpha2,scale=sigma2)))
      return(-logl)
    }
    fit = optim(initial,NegLogLike,data_obtained=data_obtained)
    alpha1r = exp(fit$par[1])
    sigma1r = exp(fit$par[2])
    alpha2r = exp(fit$par[3])
    sigma2r = exp(fit$par[4])
    pr = exp(fit$par[5])/(1 + exp(fit$par[5]))
    toreturn = c(length(fit$par),alpha1r,sigma1r,alpha2r,sigma2r,pr,-fit$value)
    return(toreturn)
  }
  else{
    NegLogLike = function(initial,data_obtained)
    {
      mu1 = initial[1]
      sigma1 = exp(initial[2])
      mu2 = initial[3]
      sigma2 = exp(initial[4])
      p = exp(initial[5])/(1 + exp(initial[5]))
      n = length(data_obtained)
      logl = sum(log(p*dlnorm(data_obtained,mu1,
                              sigma1) + (1-p)*dlnorm(data_obtained,
                                                     mean = mu2,sd=sigma2)))
      return(-logl)
    }
    fit = optim(initial,NegLogLike,data_obtained=data_obtained)
    mu1r = fit$par[1]
    sigma1r = exp(fit$par[2])
    mu2r = fit$par[3]
    sigma2r = exp(fit$par[4])
    pr = exp(fit$par[5])/(1 + exp(fit$par[5]))
    toreturn = c(length(fit$par),mu1r,sigma1r,mu2r,sigma2r,pr,-fit$value)
    return(toreturn)
  }
}

aicfunction = function(modelnumber,data_obtained)
{
  if(modelnumber == 1)
  {vals = MyMLE(data_obtained,modelnumber,
                initial = c(log(101.5127),log(1.8780),79.79213,log(6.132856),-log(1/0.35 - 1)))
  answer = 2*vals[1] - 2*vals[7]
  }
  else if(modelnumber == 2)
  {
    vals = MyMLE(data_obtained,modelnumber,initial = c(log(101.5127),log(1.8780),log(169.2757),log(2.1215),-log(1/0.35 - 1)))
    answer = 2*vals[1] - 2*vals[7]
  }
  else
  {
    vals = MyMLE(data_obtained,modelnumber,initial = c(3.9851,log(0.0098),4.3765,log(0.00589),-log(1/0.35 - 1)))
    answer = 2*vals[1] - 2*vals[7]
    
  }
  answer
}