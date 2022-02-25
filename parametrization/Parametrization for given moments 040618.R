
library(EnvStats)



SD <- 2 # Standard Deviation = sqrt(Var(X))
E <- 10 # Expected Value E(X)


# Weibull Parametrization for given E(X) and Var(x)


f <- function(k) SD^2/E^2 - (gamma(1+2/k)/gamma(1+1/k)^2)+1


parameter <- uniroot(f, c(1,1000000),tol = 0.00001)
parameter

k <- parameter$root #shape

lam <- E/gamma(1+1/k) #scale
lam

data <- rweibull(10000,k,lam)
mean(data)
sd(data)

E_wb <- lam*gamma(1+1/k)
Var_wb <- lam^2*(gamma(1+2/k)-(gamma(1+1/k))^2)
###Iteration for Weibull
f <- function(x) {
lam <- x[1]
k <- x[2]
(lam*gamma(1+1/k)-E)^2 +(lam^2*(gamma(1+2/k)-(gamma(1+1/k))^2)-SD^2)^2
}

out<-optim(c(3,3),f)
lam<- out$par[1]
k <- out$par[2]

E_wb <- lam*gamma(1+1/k)
#Check E(X) Numerically by first power Moment
f<- function(x) x*dweibull(x,scale=lam,shape=k)
integrate(f,lower=0, upper=100)

Var_wb <- lam^2*(gamma(1+2/k)-(gamma(1+1/k))^2)
###
par=c(k,lam)
pMOMENTS=function(x,para=par,pot=1){x^pot*dweibull(x,para[1],para[2])}

m1=integrate(pMOMENTS,para=par,pot=1,lower=0,upper=Inf)$value
m2=integrate(pMOMENTS,para=par,pot=2,lower=0,upper=Inf)$value

#Check Parametrization Numerically for Weibull
x=seq(1,10,,100)
mean(x)
var(x)
E<-10
Var<-4

pmoments<-function(x,para=par,pot) x^pot*para[1]/para[2]*(x/para[2])^(para[1]-1)*exp(-(x/para[2])^para[1])

funcMIN=function(para,E=E,Var=Var){
  a = para[1]; b = para[2]
  m1<- integrate(pmoments,0,Inf,para=para,pot=1)$value
  m2<- integrate(pmoments,0,Inf,para=para,pot=2)$value
  mc1=m1
  mc2=m2-mc1^2
  
  wert_E_Var=(m1-E)^2+(mc2-Var)^2
  return(wert_E_Var)
}
out<- nlminb(funcMIN,start=c(0.5,1),lower=c(0.1,0.1),upper=c(100,100),E=E,Var=Var)$par #non-finite function value

##Check moments
out
para=c(out[1],out[2])
m1<-integrate(pmoments,0,Inf,para=para,pot=1)$value #is finite
m2<-integrate(pmoments,0,Inf,para=para,pot=2)$value #is finite
mc1=m1
mc1

mc2=m2-mc1^2
mc2

##check density
f<- function(x) para[1]/para[2]*(x/para[2])^(para[1]-1)*exp(-(x/para[2])^para[1])
dweibull(2,3,1)




# Lognormal Parametrization for given E(X) and Var(X)

mu <- log(E)-0.5*log((E^2+SD^2)/E^2)
sigma2 <- log(SD^2/E^2+1, base=exp(1))
sigma <- sqrt(sigma2)

data <- rlnorm(10000,mu,sigma)
mean(data)
sd(data)

f<- function(sigma2) SD-(E*sqrt(exp(sigma2)-1))
parameter <- uniroot(f, c(0.01,10),tol = 0.00001)
sigma2<- parameter$root
sigma <- sqrt(sigma2)

mu <- log(E, base=exp(1))-0.5*sigma2
curve(f, from = -1, to = 10); abline(h = 0, lty = 3)

E_ln <- exp(mu+0.5*sigma2)
SD_ln <- E*sqrt(exp(sigma2)-1)

###Iteration for Lognormal
f <- function(x) {
  mu <- x[1]
  sigma2 <- x[2]
  (exp(mu+0.5*sigma2)-E)^2+(E*sqrt(exp(sigma2)-1)-SD)^2
}

out<-optim(c(3,3),f)
mu<- out$par[1]
sigma2 <- out$par[2]

E_ln <- exp(mu+0.5*sigma2)
SD_ln <- E*sqrt(exp(sigma2)-1)

#Check E(X) Numerically by first power Moment
f<- function(x) x*dlnorm(x,meanlog=mu,sdlog=sqrt(sigma2))
integrate(f,lower=0, upper=100)

###

# European Pareto Parametrization for given E(X) and Var(X)


f <- function(k) (E*(k-1)/k)^2*k-SD^2*(k-1)^2*(k-2)
curve(f, from = 0, to = 10)

parameter <- uniroot(f, c(0,1000000),tol = 0.00001,lower=0.0001)# lower bound must be shape > 1, else: 0/0=NA
parameter
f(6.099019)
f(1.01)
k <- parameter$root #shape

lam<- E*(k-1)/k #scale
lam
E_pareto <- k*lam/(k-1)
Var_pareto <- (E*(k-1)/k)^2*k/((k-1)^2*(k-2))

data <- rpareto(10000,lam,k)
mean(data)
sd(data)


#####
g <- function(k) (E*(k-1)/k)^2*k/((k-1)^2*(k-2))-SD^2
curve(g, from = 0, to = 10)
parameter <- uniroot(g, c(1,1000000),tol = 0.00001,lower=2) # lower bound must be >1 since the is definition gap for lam = 1 wich ends in 0/0 -> NAN instead Inf, cp 0/0=NaN and e.g. a=5 a/0=Inf
parameter
#####
###Iteration for European Pareto
f <- function(x) {
  lam <- x[1]
  k <- x[2]
  (lam*k/(k-1)-E)^2 +(lam^2*k/((k-1)^2*(k-2))-SD^2)^2
  
}

out<-optim(c(3,3),f)
lam<- out$par[1]
k <- out$par[2]

#####
###Iteration for American Pareto
E=10
SD=10.01
f <- function(x) {
  lam <- x[1]
  k <- x[2]
  (lam/(k-1)-E)^2 +(lam^2*k/((k-1)^2*(k-2))-SD^2)^2
  
}


out<-optim(c(3,3),f)
lam<- out$par[1]
k <- out$par[2]


#####


