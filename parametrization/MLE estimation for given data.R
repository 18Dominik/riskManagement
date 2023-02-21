
#data<- unlist(oprdata)

# Lognormal
###Test
data<-rlnorm(1000,2,0.2)
###

###
Fun_fit<- dlnorm
###
  ################################
  #Parameter Estimation: Maximum Likelihood by nlminb optimization
  
  loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
    par1 <- theta[1]
    par2 <- theta[2]
    
    loglik <- sum(log(Fun_fit(data, par1, par2),base=exp(1))) #FDIF
    return(-loglik)
  }
  
  
  out2 <- nlminb(c(1,1), loglik,lower=c(0.01,0.01),upper=c(1000000,1000000)) 
out2
plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2]),type="l",, main="Density Data vs. MLE")
lines(density(data),col="red")
legend("bottomright", legend=c("Data","MLE"),
       col=c("red","black"), lty=1:1,cex=0.8)
#################
  # Weibull
  ###Test
 data<-rweibull(1000,5,10)
  ###
  
  ###
  Fun_fit<- dweibull
  ###
  ################################
  #Parameter Estimation: Maximum Likelihood by nlminb optimization
  
  loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
    par1 <- theta[1]
    par2 <- theta[2]
    
    loglik <- sum(log(Fun_fit(data, par1, par2),base=exp(1))) #FDIF
    return(-loglik)
  }
  
  
  out2 <- nlminb(c(3,3), loglik,lower=c(0,0),upper=c(1000000,1000000)) 
  out2
  plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2]),type="l",, main="Density Data vs. MLE")
  lines(density(data),col="red")
  legend("bottomright", legend=c("Data","MLE"),
         col=c("red","black"), lty=1:1,cex=0.8)
  #################
  # European Stoppa II
  ###Test
 data<-rstoppa1(1000,8,6,1)
  ###
  
  ###
  Fun_fit<- dstoppa1
  ###
  ################################
  #Parameter Estimation: Maximum Likelihood by nlminb optimization
  
  loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
    par1 <- theta[1]
    par2 <- theta[2]
    par3<- theta[3]
    
    loglik <- sum(log(Fun_fit(data, par1, par2, par3),base=exp(1))) #FDIF
    return(-loglik)
  }
  
  
  out2 <- nlminb(c(1,1,1), loglik,lower=c(0.1,0.1,0.1),upper=c(1000000,1000000,1000000)) 
  out2
  plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2],out2$par[3]),type="l",, main="Density Data vs. MLE")
  lines(density(data),col="red")
  legend("bottomright", legend=c("Data","MLE"),
         col=c("red","black"), lty=1:1,cex=0.8)
  
  #################
  # European Pareto I
  ###Test
  data<-rpareto(1000,8,6)
  ###
  
  ###
  Fun_fit<- dpareto
  ###
  ################################
  #Parameter Estimation: Maximum Likelihood by nlminb optimization
  
  loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
    par1 <- theta[1]
    par2 <- theta[2]
    
    loglik <- sum(log(Fun_fit(data, par1, par2),base=exp(1))) #FDIF
    return(-loglik)
  }
  
  
  out2 <- nlminb(c(1,1), loglik,lower=c(0.001,0.001),upper=c(1000000,1000000)) 
  out2
  plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2]),type="l",, main="Density Data vs. MLE")
  lines(density(data),col="red")
  legend("bottomright", legend=c("Data","MLE"),
         col=c("red","black"), lty=1:1,cex=0.8)
  
  #################
  # Full Tails Gamma (FTG) - Castillo(2012) Full-Tails Gamma Distribution Applied to Model EV
  ###Test
 data<-rFTG(1000,50,3,14)
  ###
  
  ###
  Fun_fit<- dFTG
  ###
  ################################
  #Parameter Estimation: Maximum Likelihood by nlminb optimization
  
  loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
    par1 <- theta[1]
    par2 <- theta[2]
    par3<- theta[3]
    
    loglik <- sum(log(Fun_fit(data, par1, par2, par3),base=exp(1))) #FDIF
    return(-loglik)
  }
  
  
  out2 <- nlminb(c(1,1,1), loglik,lower=c(0.1,0.1,0.1),upper=c(1000000,1000000,1000000)) 
  out2
 plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2],out2$par[3]),type="l",, main="Density Data vs. MLE")
lines(density(data),col="red")
legend("bottomright", legend=c("Data","MLE"),
       col=c("red","black"), lty=1:1,cex=0.8)

#################

#################
# American Pareto I
###Test
data<-rpareto_2(1000,8,6)
###

###
Fun_fit<- dpareto_2
###
################################
#Parameter Estimation: Maximum Likelihood by nlminb optimization

loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
  par1 <- theta[1]
  par2 <- theta[2]
  
  loglik <- sum(log(Fun_fit(data, par1, par2),base=exp(1))) #FDIF
  return(-loglik)
}


out2 <- nlminb(c(1,1), loglik,lower=c(0.001,0.001),upper=c(1000000,1000000)) 
out2
plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2]),type="l",, main="Density Data vs. MLE")
lines(density(data),col="red")
legend("bottomright", legend=c("Data","MLE"),
       col=c("red","black"), lty=1:1,cex=0.8)

#################
# Full Tails Gamma (FTG) - Castillo(2012) Full-Tails Gamma Distribution Applied to Model EV
###Test
data<-rFTG(1000,50,3,14)
###

###
Fun_fit<- dFTG
###
################################
#Parameter Estimation: Maximum Likelihood by nlminb optimization

loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
  par1 <- theta[1]
  par2 <- theta[2]
  par3<- theta[3]
  
  loglik <- sum(log(Fun_fit(data, par1, par2, par3),base=exp(1))) #FDIF
  return(-loglik)
}


out2 <- nlminb(c(1,1,1), loglik,lower=c(0.1,0.1,0.1),upper=c(1000000,1000000,1000000)) 
out2
plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2],out2$par[3]),type="l",, main="Density Data vs. MLE")
lines(density(data),col="red")
legend("bottomright", legend=c("Data","MLE"),
       col=c("red","black"), lty=1:1,cex=0.8)

#################
# Gompertz
library(flexsurv)

###Test
data<-rgompertz(1000,8,6)
###

###
Fun_fit<- dgompertz 
###
################################
#Parameter Estimation: Maximum Likelihood by nlminb optimization

loglik <- function(theta){ #ADAPT likelihood-name and fitting parameters
  par1 <- theta[1]
  par2 <- theta[2]
  
  loglik <- sum(log(Fun_fit(data, par1, par2),base=exp(1))) #FDIF
  return(-loglik)
}


out2 <- nlminb(c(1,1), loglik,lower=c(0.001,0.001),upper=c(1000000,1000000)) 
out2
plot(sort(data),Fun_fit(sort(data),out2$par[1],out2$par[2]),type="l",, main="Density Data vs. MLE")
lines(density(data),col="red")
legend("bottomright", legend=c("Data","MLE"),
       col=c("red","black"), lty=1:1,cex=0.8)

  