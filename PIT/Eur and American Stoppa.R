
#European Stoppa

#Dichte der European Stoppa Verteilung
dstoppa1 = function(x,b,alpha, kappa){
  ifelse(x>b, kappa*alpha*(b^(alpha))*(x^(-alpha-1))*((1-(b^(alpha))*(x^(-alpha)))^(kappa-1)),0)
}

#Verteilungsfunktion der European Stoppa Verteilung
pstoppa1 = function(x,b,alpha, kappa){
  ifelse(x>b,(1-(b^(alpha))*(x^(-alpha)))^(kappa),0)
  }

#Quantilsfunktion der European Stoppa Verteilung
qstoppa1 =function(u,b,alpha, kappa){
  b*((1-u^(1/kappa))^(-1/alpha))
}

#European Stoppa verteilte Zufallszahlen (mittels Wahrscheinlichkeitsintegraltransformation)
rstoppa1=function(n,b,alpha, kappa){
  u<-runif(n)
  return(qstoppa1(u=u,b=b,alpha=alpha,kappa=kappa))
}


##
#Check PDF European Stoppa
x=seq(0,100,,200)
b=6;alpha=1;kappa=2
plot(x,dstoppa1(x,b=b, alpha=alpha,kappa=kappa),type="l",xlim=c(5,100),ylim=c(0,0.2))
lines(x,dstoppa2(x,b=b, alpha=alpha,kappa=kappa),col="blue")
legend("bottomright", legend=c("European Stoppa", "American Stoppa"), # for European Pareto, lam is a real location ("truncation") parameter
       col=c("blue", "black"), lty=1:1,cex=0.8)

integrate(dstoppa1,b=b, alpha=alpha,kappa=kappa,lower=0,upper=Inf)$value #==1?

#Check CDF European Stoppa
plot(x,pstoppa1(x,b=b, alpha=alpha,kappa=kappa),type="l",xlim=c(5,100),ylim=c(0,1))
integrate(dstoppa1,b=b, alpha=alpha,kappa=kappa,lower=b,upper=10)$value
pstoppa1(x=10,b=b,alpha=alpha,kappa=kappa) 

#Check Random Number Generation (RNG) 
testsample1=rstoppa1(10000,b=b, alpha=alpha,kappa=kappa)
plot(density(testsample1),col="red",xlim=c(0,30000))
x=seq(0,30,,200)
lines(x,dstoppa1(x,b=b,alpha=alpha,kappa=kappa),type="l",col="blue")

#Check Quantile

testvalue=20
a=pstoppa1(testvalue,b=b, alpha=alpha,kappa=kappa)
qstoppa1(a,b=b,alpha=alpha,kappa=kappa)

#Beispiel
x=10
b=5
alpha=1
kappa=2
n=4

dstoppa1(x, b, alpha, kappa)
u<-pstoppa1(x, b, alpha, kappa)
u
qstoppa1(u, b, alpha, kappa)
rstoppa1(n,b,alpha,kappa)


###
###
#American Stoppa

#Dichte der American Stoppa Verteilung
dstoppa2 = function(x,b,alpha, kappa){
  kappa*alpha*(b^(alpha))*((x+b)^(-alpha-1))*((1-(b^(alpha))*((b+x)^(-alpha)))^(kappa-1))
}

#Verteilungsfunktion der American Stoppa Verteilung
pstoppa2 = function(x,b,alpha, kappa){
  (1-(b^(alpha))*((x+b)^(-alpha)))^(kappa)
}

#Quantilsfunktion der American Stoppa Verteilung
qstoppa2 =function(u,b,alpha, kappa){
  b*((1-u^(1/kappa))^(-1/alpha))-b
}

#American Stoppa verteilte Zufallszahlen (mittels Wahrscheinlichkeitsintegraltransformation)
rstoppa2=function(n,b,alpha, kappa){
  u<-runif(n)
  return(qstoppa2(u=u,b=b,alpha=alpha,kappa=kappa))
}

#Check PDF American Stoppa
x=seq(0,100,,200)
b=6;alpha=1;kappa=2
plot(x,dstoppa1(x,b=b, alpha=alpha,kappa=kappa),type="l",xlim=c(5,100),ylim=c(0,0.2))
lines(x,dstoppa2(x,b=b, alpha=alpha,kappa=kappa),col="blue")
legend("bottomright", legend=c("European Stoppa", "American Stoppa"), # for European Pareto, lam is a real location ("truncation") parameter
       col=c("black", "blue"), lty=1:1,cex=0.8)

integrate(dstoppa2,b=b, alpha=alpha,kappa=kappa,lower=0,upper=Inf)$value #==1?

#Check CDF American Stoppa
plot(x,pstoppa2(x,b=b, alpha=alpha,kappa=kappa),type="l")
integrate(dstoppa2,b=b, alpha=alpha,kappa=kappa,lower=0,upper=20)$value
pstoppa2(x=20,b=b,alpha=alpha,kappa=kappa)

#Check Random Number Generation (RNG) 
testsample2=rstoppa2(10000,b=b, alpha=alpha,kappa=kappa)
plot(density(testsample2),col="red",xlim=c(0,30000))
x=seq(0,30,,200)
lines(x,dstoppa2(x,b=b,alpha=alpha,kappa=kappa),type="l",col="blue")

#Check Quantile

testvalue=20
a=pstoppa2(testvalue,b=b, alpha=alpha,kappa=kappa)
qstoppa2(a,b=b,alpha=alpha,kappa=kappa)

#Numeric Expectaton American Stoppa
f<- function(x) x*dstoppa2(x,b=b,alpha=alpha,kappa=kappa)
integrate(f,lower=0, upper=100000)
#Numeric Expectaton European Stoppa
f<- function(x) x*dstoppa1(x,b=b,alpha=alpha,kappa=kappa)
integrate(f,lower=0, upper=100000)

#####
#####

