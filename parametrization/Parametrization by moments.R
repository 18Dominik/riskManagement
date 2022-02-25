library(EnvStats)
E=10 # Expectation Value
Var=4 # Variance
K=250 # Kurtosis
S=10 # Skewness

#Weibull

##Given Expectation E(X) and Variance Var(X)

funcMIN=function(theta,E=E,Var=Var,K=K,S=S, t=t){
  lam = theta[1]; k = theta[2]
  m1=lam^1*gamma(1+1/k)
  m2=lam^2*gamma(1+2/k)
  m3=lam^3*gamma(1+3/k)
  m4=lam^4*gamma(1+4/k)
  mc1=m1
  mc2=m2-mc1^2
  mc3=m3-3*m2*mc1+2*mc1^3
  mc4=m4-4*m3*mc1+6*mc1^2*mc2-3*mc1^4
  
  k=mc4/(mc2)^2 #kurtosis
  s=mc3/(mc2)^1.5 #skewness
  
  wert_E_Var=(m1-E)^2+(mc2-Var)^2
  wert_K_S=(k-K)^2+(s-S)^2
  if(t==1){
  return(wert_E_Var)
    }
  if(t==2){
    return(wert_K_S)
  }
}

out<- nlminb(funcMIN,start=c(3,3),lower=c(0,0),upper=c(100,100),E=E,Var=Var,K=K,S=S,t=1)$par #t=1: Parameters for given Expectation and Variance

lam_1<- out[1] #scale
k_1<- out[2]#shape
#Check
m1=lam_1^1*gamma(1+1/k_1)
m2=lam_1^2*gamma(1+2/k_1)

mc1=m1 #=10
mc2=m2-mc1^2 #=4
mean(rweibull(10000,k_1,lam_1)) #=10
sd(rweibull(10000,k_1,lam_1))# =2 =sqrt(Var(X)=4)

out<- nlminb(funcMIN,start=c(3,3),lower=c(0,0),upper=c(100,100),E=E,Var=Var,K=K,S=S,t=2)$par#t=2 Parameters for given Kurtosis and Skewness

lam_2<- out[1] #scale
k_2<- out[2]#shape

#Check
m1=lam_2^1*gamma(1+1/k_2)
m2=lam_2^2*gamma(1+2/k_2)
m3=lam_2^3*gamma(1+3/k_2)
m4=lam_2^4*gamma(1+4/k_2)
mc1=m1
mc2=m2-mc1^2
mc3=m3-3*m2*mc1+2*mc1^3
mc4=m4-4*m3*mc1+6*mc1^2*mc2-3*mc1^4

k=mc4/(mc2)^2 #kurtosis = 250
s=mc3/(mc2)^1.5 #skewness = 10
kurtosis(rweibull(1000000,k_2,lam_2))#=250, behaves sensitive
skewness(rweibull(1000000,k_2,lam_2))


plot(density(rweibull(10000,k_1,lam_1)),xlim=c(0,100))
lines(density(rweibull(10000,k_2,lam_2)),col="red")
legend("bottomright", legend=c("Given E(X), Var(X)", "Given K(X), S(X)"),
       col=c("black", "red"), lty=1:1,cex=0.8)##Given Kurtosis K(X) and Skewness S(X)


#Lognormal

##Given Expectation E(X) and Variance Var(X)

funcMIN=function(theta,E=E,Var=Var,K=K,S=S, t=t){
  mu = theta[1]; sigma = theta[2]
  m1=exp(1*mu+1^2*sigma^2/2)
  m2=exp(2*mu+2^2*sigma^2/2)
  m3=exp(3*mu+3^2*sigma^2/2)
  m4=exp(4*mu+4^2*sigma^2/2)
  mc1=m1
  mc2=m2-mc1^2
  mc3=m3-3*m2*mc1+2*mc1^3
  mc4=m4-4*m3*mc1+6*mc1^2*mc2-3*mc1^4
  
  k=mc4/(mc2)^2 #kurtosis
  s=mc3/(mc2)^1.5 #skewness
  
  wert_E_Var=(m1-E)^2+(mc2-Var)^2
  wert_K_S=(k-K)^2+(s-S)^2
  if(t==1){
    return(wert_E_Var)
  }
  if(t==2){
    return(wert_K_S)
  }
}

out<- nlminb(funcMIN,start=c(0,0),lower=c(0,0),upper=c(100,100),E=E,Var=Var,K=K,S=S,t=1)$par #t=1: Parameters for given Expectation and Variance

mu_1<- out[1] #scale
sigma_1<- out[2]#shape

#Check
m1=exp(1*mu_1+1^2*sigma_1^2/2)
m2=exp(2*mu_1+2^2*sigma_1^2/2)
mc2=m2-m1^2
mean(rlnorm(10000,mu_1,sigma_1)) #=10
sd(rlnorm(10000,mu_1,sigma_1))#=sqrt(Var(X)=4)

# S(X) and K(X) are both only dependent only on one parameter, which is sigma i.e. sdlog
#-> cannot solve sigma for both K(X) and S(X)


plot(density(rlnorm(10000,mu_1,sigma_1)),xlim=c(0,100))



#European Stoppa II
funcMIN=function(theta,E=E,Var=Var,K=K,S=S, t=t){
  lam = theta[1]; k1 = theta[2];k2 = theta[3]
  m1=k2*beta(1-1/k1,k2)*lam^1
  m2=k2*beta(1-2/k1,k2)*lam^2
  m3=k2*beta(1-3/k1,k2)*lam^3
  m4=k2*beta(1-4/k1,k2)*lam^4

  mc1=m1
  mc2=m2-mc1^2
  mc3=m3-3*m2*mc1+2*mc1^3
  mc4=m4-4*m3*mc1+6*mc1^2*mc2-3*mc1^4
  
  k=mc4/(mc2)^2 #kurtosis
  s=mc3/(mc2)^1.5 #skewness
  
  wert_E_Var=(m1-E)^2+(mc2-Var)^2
  wert_K_S=(k-K)^2+(s-S)^2
  if(t==1){
    return(wert_E_Var)
  }
  if(t==2){
    return(wert_K_S)
  }
}

out<- nlminb(funcMIN,start=c(3,3,3),lower=c(0,0,0),upper=c(100,100,100),E=E,Var=Var,K=K,S=S,t=1)$par #t=1: Parameters for given Expectation and Variance

lam_1<- out[1] #scale
k1_1<- out[2]#shape1
k2_1<-out[3]#shape2

#Check:
m1=k2_1*beta(1-1/k1_1,k2_1)*lam_1^1#=10
m2=k2_1*beta(1-2/k1_1,k2_1)*lam_1^2
mc2=m2-m1^2 #=4
set.seed(1)
mean(rstoppa1(1000000,lam_1,k1_1,k2_1))#=10
sd(rstoppa1(1000000,lam_1,k1_1,k2_1))#=2 =sqrt(Var(X)=4)


mean(rstoppa2(1000000,lam_1,k1_1,k2_1))#=10
sd(rstoppa2(1000000,lam_1,k1_1,k2_1))#=2 =sqrt(Var(X)=4)
lam_1

data_eu<-rstoppa1(1000000,lam_1,k1_1,k2_1)
data_ae<-rstoppa2(1000000,lam_1,k1_1,k2_1)
plot(density(data_eu),xlim=c(-1,20))
lines(density(data_ae),col="red")

out<- nlminb(funcMIN,start=c(0,0,0),lower=c(0,0,0),upper=c(100,100,100),E=E,Var=Var,K=K,S=S,t=2)$par #t=1: Parameters for given Expectation and Variance




#European Pareto
funcMIN=function(theta,E=E,Var=Var,K=K,S=S, t=t){
  lam = theta[1]; k = theta[2]
  m1=k/(k-1)*lam^1
  m2=k/(k-2)*lam^2
  m3=k/(k-3)*lam^3
  m4=k/(k-4)*lam^4
  
  mc1=m1
  mc2=m2-mc1^2
  mc3=m3-3*m2*mc1+2*mc1^3
  mc4=m4-4*m3*mc1+6*mc1^2*mc2-3*mc1^4
  
  k=mc4/(mc2)^2 #kurtosis
  s=mc3/(mc2)^1.5 #skewness
  
  wert_E_Var=(m1-E)^2+(mc2-Var)^2
  wert_K_S=(k-K)^2+(s-S)^2
  if(t==1){
    return(wert_E_Var)
  }
  if(t==2){
    return(wert_K_S)
  }
}

out<- nlminb(funcMIN,start=c(3,3),lower=c(0,0),upper=c(100,100),E=E,Var=Var,K=K,S=S,t=1)$par #t=1: Parameters for given Expectation and Variance
lam_1<- out[1]
k_1<- out[2]

mean(rpareto(10000,lam_1,k_1))#=10
sd(rpareto(10000,lam_1,k_1))#=2 =sqrt(Var(X)=4)
mean(rpareto_2(10000,lam_1,k_1))
lam_1
k_1

data_eu<-rpareto(1000000,lam_1,k_1)
data_ae<-rpareto2(1000000,lam_1,k_1)
plot(density(data_eu),xlim=c(-1,20))
lines(density(data_ae),col="red")


