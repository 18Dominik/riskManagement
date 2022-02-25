#alpha = 0.05
#DATA = rt(2, 50)
#DATA = DATA[order(DATA)]
#base_mean = mean(DATA)
#base_sd= sd(DATA)
#R=10
#TN_B(data, alpha, mean(data), sd(data), pnorm, T=0.0,R=100)

DEBUG <- FALSE
myprint <- if(DEBUG) print else function(x){}

TN_B =function(DATA=data, p=alpha, base_mean=mean(data), base_sd=sd(data), FUN="pnorm", distn="pnorm", T=0.0,R=100){
  ld=length(DATA)

  hist(DATA)

  ECDF=(1:ld)/ld                                                                                                # Calculate empirical cumulative distribution function
  FCDF= pnorm(DATA, mean=base_mean, sd=base_sd)                    # Import fitted values in loss space
  myprint("FCDF"); myprint(head(FCDF))


  plot(DATA, ECDF, type="l", col="blue")                               # LOSS SPACE
  points(DATA, FCDF, type = "l",col="red")

  plot(ECDF, ECDF, type = "l",col="blue")                               # PROBABILITY SPACE
  points(FCDF, ECDF, type = "l", col="red")

  W= sqrt(2)*(ECDF[2:ld]-ECDF[1:ld-1]) #distance between adjacent points on 45°line is always the same
  myprint("W"); myprint(head(W))
  H= abs(FCDF-ECDF)/sqrt(2)
  lH=length(H)
  plot(1:ld,H, type="l")
  A=abs((H[1:lH-1]+H[2:lH]/2)*W[1:lH-1])  # Trapezium Areas
  A_total=sum(A)
  myprint("A_total"); myprint(A_total)


  # RESAMPLING PART

  set_all = cbind(base_mean*seq(0, 2, 2/R),base_sd*seq(0, 2, 2/R)) # Dominik: Sequenz geändert
  resample_par=cbind(sample(set_all[,1],1+R, replace=TRUE),sample(set_all[,2],1+R, replace=TRUE))
  myprint("resample_par"); myprint(resample_par)


  H = FCDFS = matrix(, nrow=1+R, ncol=ld)
  for (k in 1: (1+R)){
    FCDFS[k,]= FUN(DATA, mean=resample_par[k, 1], sd=resample_par[k,2]) # Dominik: hier DATA statt ECDF
    H[k,]= abs(FCDFS[k,]-ECDF)/sqrt(2)
  }
  myprint("H"); myprint(head(H[,1:5]))

  A = matrix(, nrow = R+1, ncol = lH-1)
  for (k in 1:(R+1)){
    A[k,]=abs((H[k,1:(lH-1)]+H[k,2:lH]/2)*W[1:(lH-1)])
  }
  myprint("A"); myprint(head(A[,1:5]))

  A_k=apply(A,1,sum)                     # Summation über Spalten  CHECK !!!!!
  A_mean=mean(A_k)
  myprint("A_k"); myprint(head(A_k))

  myprint("A_mean"); myprint(A_mean)
  difference = abs(A_k-A_mean)   # Dominik: Betrag der Differenz abs() !!!
  myprint("difference"); myprint(head(difference))

  count_difference=length(which(difference>A_total))
  myprint("count_difference"); myprint(head(count_difference))
  # calculate p-value
  p_TNB <- (1+count_difference)/(2+R)
  if (p_TNB > p){Message="Accept H_0"}
  else {Message="Reject H_0"}

  return(list(Decision=Message,p=p_TNB))
}

# Test function

  alpha = 0.05
  set.seed(1)

  data = rt(5000, 20)
  data = data[order(data)]
  load("TNBTestdata.Rdata")
  TN_B(data, alpha, mean(data), sd(data), pnorm, T=0.0,R=100)

#system.time(source("TN_B_Steffan_Fischer_final.R"))
