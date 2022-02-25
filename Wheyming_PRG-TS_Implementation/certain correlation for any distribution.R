#based on paper Wheyming (1996):Generating pseudo-random time series with specifiedmarginal distributions
n=10000
m=10#number of iterations
data<-numeric(n)
corrx<- numeric(m)
corry<- numeric(m)
targety<-0.714 #set target corr
rx<-targety #Step 0
data_lag<-numeric(n)

var_res<-1-rx^2
sigma_res<-sqrt(var_res)
pFUN<-pweibull
qFUN<-qweibull
par1<-1
par2<-1

shift<- function(data){
  for(i in 2:n){
    data_lag[i]<- data[i-1]
    
  }
  return (data_lag)
}

#Step 1
auto<-function(){
  old.rx<-targety
repeat{
  
for(j in 1:m){
  data[1]<-0
  for(i in 2:n){
    data[i]<- old.rx*data[i-1]+rnorm(1,0,sigma_res)
    
  }

 
  data_lag<-shift(data)
  corrx[j]<- cor(data,data_lag) # = rx = given correlation
  
  u<-pFUN(data,par1,par2)
  y<-qFUN(u,par1,par2)
  
  y_lag<-shift(y)
  corry[j]<- cor(y,y_lag) 
}
plot(density(data))
#Step 2
corry_mean=sum(corry)/m
se_corry_mean<- sqrt(m^(-1)*(m-1)^(-1)*sum((corry-corry_mean)^2))
#Step 3

if(abs(corry_mean-targety)>se_corry_mean){
  delta<- abs(corry_mean-targety)
  if(corry_mean<targety){
    new.rx<-old.rx+delta
  }else if(corry_mean>targety){
    new.rx<-old.rx-delta}
  print(corry_mean-targety)
}else if(abs(corry_mean-targety)<=se_corry_mean) return(old.rx)
old.rx<- new.rx
}
}
corr_it<-auto()

#check
data[1]<-0
for(i in 2:n){
  data[i]<- corr_it*data[i-1]+rnorm(1,0,sigma_res)
  
}
data_lag<-shift(data)
corrx<- cor(data,data_lag) # = rx = given correlation

u<-pFUN(data,par1,par2)
y<-qFUN(u,par1,par2)

y_lag<-shift(y)
corry<- cor(y,y_lag) 
corry #=target corr
