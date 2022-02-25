TN_A =function(DATA=data, p=alpha, base_mean=mean(data), base_sd =sd(data), FUN=pnorm, T=0.0){

  # base_sd and base_mean are the parameters of the fitting curve
  # T = Threshold for left truncated loss data
  DATA <-tail(DATA,(1-T)*length(DATA)) # Left truncation if T>0
  ld=length(DATA)

  hist(DATA)

  ECDF=(1:ld)/ld                                                                                                 # Calculate empirical cumulative distribution function
  FCDF= FUN(DATA, mean=base_mean, sd=base_sd)                         # Import fitted values in loss space

  # Plot in Loss Space

  plot(DATA, ECDF, type="l", col="blue")
  points(DATA, FCDF, type = "l",col="red")

  # Plot in Probability Space

  plot(ECDF, ECDF, type = "l",col="blue")
  points(FCDF, ECDF, type = "l", col="red")

  # Calculate W[i] as distance of Y[i+1] and [Y_i] on the 45Â° line

  W= sqrt(2)*(ECDF[2:ld]-ECDF[1:ld-1]) #distance between adjacent points on 45Â°line is always the same

  # Calculate H as distance between transformed fitting distribution and 45Â° line

  H= abs(FCDF-ECDF)/sqrt(2)
  lH=length(H)

  plot(1:ld,H, type="l")

  # Calculate area between transoformed fitting distribtuion and 45Â° line by calcualting trapzezium area (a+c)/2 *h ) m*h = A

  A=abs((H[1:lH-1]+H[2:lH]/2)*W[1:lH-1])  # Trapezium Area

  # Calculate total discrepance between tranformed fitting distribution and 45Â° line
  A_total=sum(A)

  # TN-A-test for enclosed Area with H_0: A_total = 0
  ## calculate critical area A_c  for significance level p

  A_c = 2*sqrt(2)*p*(1-sqrt(2)*p-T) # cp. Mitic (2015), p. 97

  if (A_total<A_c){Message="Accept H_0"}
  else {Message="Reject H_0"}

  return(list(Decision=Message,AreaData=A_total,CriticalArea=A_c))

}



alpha = 0.05
data = rt(5000, 50)
data = data[order(data)]
TN_A(data, alpha, mean(data), sd(data), pnorm, T=0.0)

