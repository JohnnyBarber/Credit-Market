###################################################################
#                          Problem 9                              #
###################################################################
library(cubature)

#-----------------------------b-----------------------------------#
f_zero = function(K,theta,sig,lambda0,r,t){
  gamma = sqrt(K^2+2*sig^2)
  lead = 2*K*theta/(sig^2)
  A = lead*log(2*gamma*exp((K+gamma)*t/2)/(2*gamma+(K+gamma)*(exp(gamma*t)-1)))
  B = 2*(exp(gamma*t)-1)/(2*gamma+(K+gamma)*(exp(gamma*t)-1))
  result = exp(A-B*lambda0)*exp(-r*t)
  return(result)
}

K = 1; theta = 0.02; sig = 0.15; lambda0 = 0.01; r = 0.01
bonds = seq(10)
yields = seq(10)
spreads = seq(10)

for(i in 1:10){
  bonds[i] = f_zero(K,theta,sig,lambda0,r,i)
  yields[i] = -log(bonds[i])/i
  spreads[i] = yields[i] - r
}

#-----------------------------c-----------------------------------#
limit_bond = f_zero(K,theta,sig,lambda0,r,0.0000001)
limit_yield = -log(limit_bond)/0.0000001
limit_spread = limit_yield - r

#-----------------------------d-----------------------------------#
freq = 0.5
R = 0.5
c = 0.02

f_coupon = function(K,theta,sig,lambda0,r,t,freq,c){
  result = 0
  period = t/freq
  for(i in 1:period){
    result = result + exp(-r*K*freq)*c*freq*f_zero(K,theta,sig,lambda0,r,t)
  }
  return(result)
}

f_recovery = function(K,theta,sig,lambda0,r,t,R){
  f_temp = function(t){
    result = (-r)*exp(-r*t)*(1-f_zero(K,theta,sig,lambda0,r,t))
    return(result)
  }
  integral = adaptIntegrate(f_temp,lower = 0, upper = t)$integral
  recovery = R*(exp(-r*t)*(1-f_zero(K,theta,sig,lambda0,r,t))-integral)
  return(recovery)
}

f_notional = function(K,theta,sig,lambda0,r,t){
  return(exp(-r*t)*f_zero(K,theta,sig,lambda0,r,t))
}

price = seq(10)
for(i in 1:10){
  price[i] = f_coupon(K,theta,sig,lambda0,r,i,freq,c) + 
    f_recovery(K,theta,sig,lambda0,r,i,R) + 
    f_notional(K,theta,sig,lambda0,r,i)
}


