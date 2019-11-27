####Credit Risk HW4####
#-------------------P10---------------------#
kappa=0.05
theta=0.03
sigma=0.1
lambda0=0.02
gamma=sqrt(kappa^2+2*sigma^2)
r=0.01
R=0.5

AB<-function(t){
  a1=2*kappa*theta/(sigma^2)
  a2=2*gamma*exp((kappa+gamma)*(t/2))
  a3=2*gamma+(kappa+gamma)*(exp(gamma*t)-1)
  A=a1*log(a2/a3)
  b1=2*(exp(gamma*t)-1)
  b2=2*gamma+(kappa+gamma)*(exp(gamma*t)-1)
  B=b1/b2
  CIR=exp(A-B*lambda0)
  return(CIR)
}
dis<-function(t){
  disc<-exp(-r*t)
  return(disc)
}

CH<-function(t){
  c1=-(2*kappa*theta)/(sigma^2)*(kappa^2-gamma^2)*(exp(gamma*t)-1)
  c2=4*gamma+2*(kappa+gamma)*(exp(gamma*t)-1)
  C=c1/c2
  h1=4*gamma^2*exp(gamma*t)
  h2=(2*gamma+(kappa+gamma)*(exp(gamma*t)-1))^2
  H=h1/h2
  out=C+H*lambda0
  return(out)
}

#-------------(a) & (b)---------------#
library(cubature)

part = function(s){
  return(dis(s)*CH(s)*AB(s))
}

#for discrete: quarterlay fee payments 
sprd_d=c()
for (s in 1:10){
  t_d=seq(0.25,s,by=0.25)
  V0fee_d=sum(0.25*dis(t_d)*AB(t_d))
  V0pro_d=(1-R)*adaptIntegrate(part,lowerLimit = 0, upperLimit = s)$integral
  sprd_d[s]=V0pro_d/V0fee_d
}

#--------------(c)-------------------#
r=0.05
sprd_n=c()

for (s in 1:10){
  t_d=seq(0.25,s,by=0.25)
  V0fee_n=sum(0.25*dis(t_d)*AB(t_d))
  V0pro_n=(1-R)*adaptIntegrate(part,lowerLimit = 0, upperLimit = s)$integral
  sprd_n[s]=V0pro_n/V0fee_n
}

#--------------(d)-------------------#
#for continous
first_part = function(s){
  return(dis(s)*AB(s))
}

sprd_c=c()
for (s in 1:10){
  V0fee_c=adaptIntegrate(first_part,lowerLimit = 0, upperLimit = s)$integral
  V0pro_c=(1-R)*adaptIntegrate(part,lowerLimit = 0, upperLimit = s)$integral
  sprd_c[s]=V0pro_c/V0fee_c
}
