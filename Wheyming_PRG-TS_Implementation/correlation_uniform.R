#based on paper R.Thomas Willemain & A.Philip Desautels: A method to generate autocorrelated uniform random numbers,1993
#Journal of statistical computation and simulation, 45:1-2,23-31, DOI: 10.1080/00949659308811469

c=0.5
n=100
x=numeric(n)
x_lag1= numeric(n-1)


psum_unif <- function(x){
  if(c>=1){
    if(x>=0 && x<=1) {psum=0.5*x^2/c
    }else if(x>1 && x<=c){psum= (x-0.5)/c
    }else if (x>c && x<= c+1){psum= 1-0.5*(1+c-x)^2/c
    }else{print("NA c>=1")}
 
  }else if(c<1 && c>=0){
    if(x>=0 && x <=c){psum=0.5*x^2/c
    }else if(x>c && x<=1){psum=x-0.5*c
    }else if (x>1 && x<=c+1){psum=(1-0.5*c)+(1+1/c)*(x-1)-0.5*(x^2-1)/c
    }else{print("NA c<1")}
  }
  return (psum)
  }

auto<- function(corr=corr){
  x[1]<- runif(1)+c*runif(1)# random number U(0,1)
  
  if(corr==1){ #positive correlation
    for(i in 2:n){
      x[i]<- psum_unif(x[i-1])+c*runif(1)
    }
    }else if(corr==0){ #negative correlation
      for(i in 2:n){
        x[i]<- 1-psum_unif(x[i-1])+c*runif(1)
      }
}
return (x)
}

shift<- function(x){
  for(i in 2:n){
    x_lag1[i]<- x[i-1]
  }
  return (x_lag1)
}
corr0=0 #negative correlation
corr1=1 #positive correlation
x0 <-auto(corr=corr0)
x1<-auto(corr1)

#check
max(x0 )
max(x1)
min(x1)


plot(x0 ,ylim=c(-0.1,0.7),xlim=c(1,n),type="l")
lines(x1,type="l",col="red")
u1<- punif(x0 )
u2<- punif(x1)
plot(u1,ylim=c(0,1),xlim=c(1,n),type="l")
lines(u2,type="l",col="red")

x_norm1<- dnorm(x0 ,0,1)
x_norm2<- dnorm(x1,0,1)

plot(density(x_norm1))
lines(density(x_norm2),col="red")
plot(x_norm1,ylim=c(0,1),xlim=c(1,n),type="l")
lines(x_norm2,type="l",col="red")

plot(density(x0 ))
#check for c=100 -> lim_corr(c->inf)=0
plot(x0 )
plot(x1)

#positive correlation
x1_lag<- shift(x1)
plot(x1_lag,x1,type="p")
plot(x1_lag)

#negative correlation
x0_lag<- shift(x0 )
plot(x0_lag,x0 ,type="p")
plot(x0_lag)

cor(x0,x0_lag)
cor(x1,x1_lag)

