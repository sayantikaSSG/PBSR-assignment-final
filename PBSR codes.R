#3
library(scales)
attach(faithful)
hist(faithful$waiting,xlab = 'waiting',probability = T,col='pink',main='')
data_q3 <- faithful 
x <- sort(data_q3$waiting)

p <- length(x[x<65])/length(x)
as <- mean(x[x<65])
ass <- var(x[x<65])
s <- ass/as
a <- as/s
mu <- mean(x[x>=65])
sigma <- sd(x[x>=65])
theta_inital <- c(p, a, s, mu, sigma)
neg_log_likelihood <- function(theta, data){
  n = length(data)
  
  p = theta[1]
  a = theta[2]
  s = theta[3]
  mu = theta[4]
  sigma = theta[5]
  
  l = 0
  for (i in 1:n) {
    l = l + log(p*dgamma(data[i], shape = a, scale = s) + (1-p)*dnorm(data[i], mean = mu, sd = sigma))
  }
  return(-l)
}
fit = optim(theta_inital,
            neg_log_likelihood,
            data = x,
            control = list(maxit = 1500),
            lower = c(0, 0, 0, -Inf, 0),
            upper = c(1, Inf, Inf, Inf, Inf),
            method="L-BFGS-B")
theta_1 = fit$par
theta_1
p = theta_1[1]
a = theta_1[2]
s = theta_1[3]
mu = theta_1[4]
sigma = theta_1[5]
model_1 = p*dgamma(x, shape = a, scale = s) + (1-p)*dnorm(x, mean = mu, sd = sigma)
aic_1 <- 2*5 + neg_log_likelihood(theta_1, x)

# model 2
p <- length(x[x<65])/length(x)
as_1 <- mean(x[x<65])
ass_1 <- var(x[x<65])
s_1 <- ass_1/as_1
a_1 <- as_1/s_1
as_2 <- mean(x[x>=65])
ass_2 <- var(x[x>=65])
s_2 <- ass_2/as_2
a_2 <- as_2/s_2
theta_inital <- c(p, a_1, s_1, a_2, s_2)
neg_log_likelihood <- function(theta, data){
  n <- length(data)
  
  p <- theta[1]
  a_1 <- theta[2]
  s_1 <- theta[3]
  a_2 <- theta[4]
  s_2 <- theta[5]
  
  l <- 0
  for (i in 1:n) {
    l = l + log(p*dgamma(data[i], shape = a_1, scale = s_1) + (1-p)*dgamma(data[i], shape = a_2, scale = s_2))
  }
  return(-l)
}
fit = optim(theta_inital,
            neg_log_likelihood,
            data = x,
            control = list(maxit = 1500),
            lower = c(0, 0, 0, 0, 0),
            upper = c(1, Inf, Inf, Inf, Inf),
            method="L-BFGS-B")
theta_2 <- fit$par
theta_2
p <- theta_2[1]
a_1 <- theta_2[2]
s_1 <- theta_2[3]
a_2 <- theta_2[4]
s_2 <- theta_2[5]
model_2 <- p*dgamma(x, shape = a_1, scale = s_1) + (1-p)*dgamma(x, shape = a_2, scale = s_2)
aic_2 <- 2*5 + neg_log_likelihood(theta_2, x)

# model 3
p <- length(x[x<65])/length(x)
m_1 <- mean(x[x<65])
v_1 <- var(x[x<65])
sigma2_1 <- log((v_1/m_1^2) + 1)
mu_1 <- log(m_1) - sigma2_1/2
m_2 <- mean(x[x>=65])
v_2 <- var(x[x>=65])
sigma2_2 <- log((v_2/m_2^2) + 1)
mu_2 <- log(m_2) - sigma2_2/2
theta_inital <- c(p, mu_1, sqrt(sigma2_1), mu_2, sqrt(sigma2_2))
neg_log_likelihood <- function(theta, data) {
  n <- length(data)
  
  p <- theta[1]
  mu_1 <- theta[2]
  sigma_1 <- theta[3]
  mu_2 <- theta[4]
  sigma_2 <- theta[5]
  
  l <- 0
  for (i in 1:n) {
    l = l + log(p*dlnorm(data[i], meanlog = mu_1, sdlog = sigma_1) + (1-p)*dlnorm(data[i], meanlog = mu_2, sdlog = sigma_2))
  }
  
  return(-l)
}
fit = optim(theta_inital,
            neg_log_likelihood,
            data = x,
            control = list(maxit = 1500),
            lower = c(0, -Inf, 0, -Inf, 0),
            upper = c(1, Inf, Inf, Inf, Inf),
            method="L-BFGS-B")
theta_3 <- fit$par
theta_3
p <- theta_3[1]
mu_1 <- theta_3[2]
sigma_1 <- theta_3[3]
mu_2 <- theta_3[4]
sigma_2 <- theta_3[5]
model_3 <- p*dlnorm(x, meanlog = mu_1, sdlog = sigma_1) + (1-p)*dlnorm(x, meanlog = mu_2, sdlog = sigma_2)
aic_3 <- 2*5 + neg_log_likelihood(theta_3, x)

hist(x, xlab = 'waiting', probability = T, col='pink', main='')
lines(x, model_1, col = "red")
lines(x, model_2, col = "blue")
lines(x, model_3, col = "green")

results <- data.frame(
  models = c("Gamma + Normal", "Gamma + Gamma", "Lognormal + Lognormal"),
  AIC = c(aic_1, aic_2, aic_3)
)
results

p <- theta_3[1]
mu_1 <- theta_3[2]
sigma_1 <- theta_3[3]
mu_2 <- theta_3[4]
sigma_2 <- theta_3[5]
reqd_prob <- (p*plnorm(70, meanlog = mu_1, sdlog = sigma_1) + (1-p)*plnorm(70, meanlog = mu_2, sdlog = sigma_2)) - (p*plnorm(60, meanlog = mu_1, sdlog = sigma_1) + (1-p)*plnorm(60, meanlog = mu_2, sdlog = sigma_2))
reqd_prob

##==============================================================================================================
##==============================================================================================================
##Problem 4
library(ggplot2)
library(MASS)
library(stats)
plot(Insurance$Holders,Insurance$Claims
     ,xlab = 'Holders',ylab='Claims',pch=20)
grid()


############   PART A   ####################
#Claimsi = β0 + β1 Holdersi + εi, i = 1, 2, · · · , n

# Claimsi ∼ N(μi, σ2), where
#μi = β0 + β1 Holdersi + εi, i = 1, 2, · · · , n


# Neg log likelihood

NegLogLike<- function(theta,data){
  sigma = theta[3]
  n = nrow(data)
  l=0
  for(i in 1:n){
    mu = theta[1]+theta[2]*data$Holders[i]
    l = l + log(dnorm(data$Claims[i],mean = mu,sd=sigma))
  }
  return(-l)
}


theta_initial=c(0.01,0.1,10)
NegLogLike(theta_initial,Insurance)

fit = optim(theta_initial
            ,NegLogLike
            ,data=Insurance,
            lower = 0, upper = Inf)
ggplot(data=Insurance)+geom_line(aes(Holders, fit$par[1]+fit$par[2]*Holders))+geom_point(aes(Holders,Claims))
theta_hat = fit$par
theta_hat



BIC=2*NegLogLike(theta_hat,Insurance)+log(nrow(Insurance))*2

BIC


#---------------------------------------------------------#

#### PART B ######


#Claimsi = β0 + β1 Holdersi + εi, i = 1, 2, · · · , n
#Assume : εi ∼ Laplace(0, σ2). Note that β0, β1 ∈ R and σ ∈ R+.
#(i) Clearly write down the negative-log-likelihood function in R. Then use optim function to estimate MLE
##of θ = (β0, β1, σ)
#(ii) Calculate Bayesian Information Criterion (BIC) for the model.


laplace<-function(x,loc, b){
  exp(-abs(x-loc)/b)/(2*b)
}


NegLogLike<- function(theta,data){
  sigma = theta[3]
  n = nrow(data)
  l=0
  for(i in 1:n){
    mu = theta[1]+theta[2]*data$Holders[i]
    l = l + log(laplace(data$Claims[i],loc = mu,b = sigma))
    #print(l)
  }
  return(-l)
}



theta_initial=c(8,0.1,10)
NegLogLike(theta_initial,Insurance)

fit = optim(theta_initial
            ,NegLogLike
            ,data=Insurance)


ggplot(data=Insurance)+geom_line(aes(Holders, fit$par[1]+fit$par[2]*Holders))+geom_point(aes(Holders,Claims))
theta_hat = fit$par
theta_hat



BIC=2*NegLogLike(theta_hat,Insurance)+log(nrow(Insurance))*2
BIC




####PART C



dataIn=Insurance[-c(61),]
NegLogLike<- function(theta,data){
  sigma = theta[3]
  n = nrow(data)
  l=0
  for(i in 1:n){
    mu = theta[1]+theta[2]*data$Holders[i]
    l = l + log(dlnorm(data$Claims[i],meanlog = mu,sdlog = sigma))
    
  }
  return(-l)
}


theta_initial=c(1,0.1,10)
NegLogLike(theta_initial,dataIn)

fit = optim(theta_initial
            ,NegLogLike
            ,data=dataIn)

ggplot(data=dataIn)+geom_line(aes(Holders, fit$par[1]+fit$par[2]*Holders))+
  geom_point(aes(Holders,Claims))
theta_hat = fit$par
theta_hat



BIC=2*NegLogLike(theta_hat,dataIn)+log(nrow(dataIn))*2
BIC


## PART D



NegLogLike<- function(theta,data){
  sigma = theta[3]
  n = nrow(data)
  l=0
  for(i in 1:n){
    if (data$Claims[i]!=0){
    mu = theta[1]+theta[2]*data$Holders[i]
    l = l + log(dgamma(data$Claims[i],shape = mu,scale = sigma))
    #print(l)
    
    }
  }
  return(-l)
}

theta_initial=c(1,0.1,1)
NegLogLike(theta_initial,Insurance)


fit = optim(theta_initial
            ,NegLogLike
            ,data=Insurance)


theta_hat = fit$par
theta_hat
BIC=2*NegLogLike(theta_hat,Insurance)+log(nrow(Insurance))*2
BIC


#=======================================================================================================================================================


##==============================================================================================================
##==============================================================================================================
#5
library(quantmod)
getSymbols('TCS.NS')
tail(TCS.NS)
plot(TCS.NS$TCS.NS.Adjusted)
getSymbols('^NSEI')
tail(NSEI)
plot(NSEI$NSEI.Adjusted)
TCS_rt = diff(log(TCS.NS$TCS.NS.Adjusted))
Nifty_rt = diff(log(NSEI$NSEI.Adjusted))
retrn = cbind.xts(TCS_rt,Nifty_rt)
retrn = na.omit(data.frame(retrn))
plot(retrn$NSEI.Adjusted,retrn$TCS.NS.Adjusted
     ,pch=20
     ,xlab='Market Return'
     ,ylab='TCS Return'
     ,xlim=c(-0.18,0.18)
     ,ylim=c(-0.18,0.18))
grid(col='grey',lty=1)
#rtcs = a + b rni + error
rtcs = mean(retrn$TCS.NS.Adjusted)
rni = mean(retrn$NSEI.Adjusted)
row = cor(retrn$TCS.NS.Adjusted,retrn$NSEI.Adjusted)
sdtcs = sqrt(var(retrn$TCS.NS.Adjusted))
sdnifty = sqrt(var(retrn$NSEI.Adjusted))
n = nrow(retrn)
b1 = row*(sdtcs/sdnifty)
a1 = sdtcs - b1*sdnifty
lin_mod = summary(lm(TCS.NS.Adjusted~NSEI.Adjusted, data = retrn))
a2 = lin_mod$coefficients [1,1]
b2 = lin_mod$coefficients [2,1]
rtcs_hat = a1 + b1 * retrn$NSEI.Adjusted
error = retrn$TCS.NS.Adjusted - rtcs_hat
sigma1 = sqrt(var(error))
Method_of_Moments <- c(a1,b1,sigma1)
k = (n-2)/n
errornew = retrn$TCS.NS.Adjusted - (a2+ b2*retrn$NSEI.Adjusted)
sigma2 = (sqrt(var(errornew)))*k
OLS <- c(a2,b2, sigma2)
Parameters <- c("alpha", "beta", "sigma")
table = data.frame(Parameters, Method_of_Moments,OLS)
table
est1 = 3200 - ( 200*b1)
est2 = 3200 - ( 200*b2)


