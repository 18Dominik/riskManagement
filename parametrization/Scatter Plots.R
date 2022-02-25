par(mfrow=c(1,2))

#Scatter Plot E(X)~Var(X) European Pareto
set.seed(1)
k<-runif(1000000,2,100)
lam<-runif(1000000,0,100)
#lam<- 8.36
#k<- 6.099
E<-lam*k/(k-1)
Var<-(lam^2*k)/((k-1)^2*(k-2))
plot(E,Var,ylim=c(0,10),xlim=c(0,40),main="European Pareto I, shape in 2:100 and location in 0:100",xlab="Expectation Value",ylab="Variance") 

#Scatter Plot E(X)~Var(X) American Pareto
# Does not work for given paramter E(X)=10 and SD(X)=2. SD(X) must be >10
#for lim SD(X) -> 10, lim k -> +Inf
set.seed(1)
k<-runif(1000000,2,100)
lam<-runif(1000000,0,100)
#lam<- 8.3603
#k<- 6.099
E<-lam/(k-1)
Var<-(lam^2*k)/((k-1)^2*(k-2))
plot(E,Var,ylim=c(0,100),xlim=c(0,40),main="American Pareto I, shape in 2:100 and location in 0:100",xlab="Expectation Value",ylab="Variance") 




#Scatter Plot E(X)~Var(X) Lognormal
set.seed(1)
mu<-runif(1000000,-100,100)
sigma<-runif(1000000,0,100)
#mu<- 2.2830
#sigma<- 0.198
E<-exp(mu+sigma^2/2)
Var<-(exp(sigma^2)-1)*exp(2*mu+sigma^2)
plot(E,Var,ylim=c(0,100),xlim=c(0,40),main="Scatter-Plot for Lognormal Distribution for mu in -100:100 and sigma in 0:100",xlab="Expectation Value",ylab="Variance") 

#Scatter Plot E(X)~Var(X) Weibull
set.seed(1)
lam<-runif(10000,0,100)
k<-runif(10000,0,100)
#lam<- 10.7998
#k<- 5.7974
E<-lam*gamma(1+1/k)
Var<-lam^2*(gamma(1+2/k)-(gamma(1+1/k))^2)
plot(E,Var,ylim=c(0,100),xlim=c(0,40),main="Scatter-Plot for Weibull Distribution for shape and scale in 0:100",xlab="Expectation Value",ylab="Variance") 


#Scatter Plot E(X)~Var(X) European Stoppa
set.seed(1)
lam<-runif(1000000,0,100)
k1<-runif(1000000,0,100)
k2<- runif(1000000,0,100)
#lam<- 8.424185
#k1<- 6.018343
#k2<- 0.9121991


m1=k2*beta(1-1/k1,k2)*lam^1
m2=k2*beta(1-2/k1,k2)*lam^2
E<-m1
Var<-m2-m1^2


plot(E,Var,ylim=c(0,100),xlim=c(0,40),main="Scatter-Plot for European Stoppa II Distribution for shape1/shape2 and scale in 0:100",xlab="Expectation Value",ylab="Variance") 

#Scatter Plot E(X)~Var(X) FTG
set.seed(1)
rho<-runif(100000,0,100)
alpha<-runif(100000,-100,100)
theta<- runif(100000,0,100)
#rho<- 14.64466
#alpha<- 50
#theta<- 3.535534
library(expint)

mu<- exp(-rho)*rho^(alpha)/gammainc(alpha,rho )
E<- (alpha-rho+mu)/theta
Var<- (alpha+(1+rho-alpha)*mu-mu^2)/theta^2


plot(E,Var,ylim=c(0,100),xlim=c(0,40),main="Scatter-Plot for FTG Distribution for alpha in -100:100, rho and theta in 0:100",xlab="Expectation Value",ylab="Variance") 



