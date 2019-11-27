####Credit Risk HW4####
rm(list=ls())
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

compute_vport = function(t){
  return(dis(t)*CH(t)*AB(t))
}

#for discrete: quarterlay fee payments 
sprd_d=c()
for (s in 1:10){
  t_d=seq(0.25,s,by=0.25)
  V0fee_d=sum(0.25*dis(t_d)*AB(t_d))
  V0pro_d=(1-R)*adaptIntegrate(compute_vport,lowerLimit = 0, upperLimit = s)$integral
  sprd_d[s]=V0pro_d/V0fee_d
}

#--------------(c)-------------------#
r_new=0.05
sprd_n=c()
for (T in 1:10){
  t_d=seq(0.25,T,by=0.25)
  V0fee_n=sum(dis(t_d,r_new)*AB(t_d,lambda0))
  V0pro_n=(1-R)*sum(dis(t_d,r_new)*CH(t_d,lambda0)*AB(t_d,lambda0))
  sprd_n[T]=V0pro_n/V0fee_n
}

#--------------(d)-------------------#
#for continous
delta=1/10000
sprd_c=c()
for (T in 1:10){
  t_c=seq(delta,T,by=delta)
  V0fee_c=sum(dis(t_c,r)*AB(t_c,lambda0))
  V0pro_c=(1-R)*sum(dis(t_c,r)*CH(t_c,lambda0)*AB(t_c,lambda0))
  sprd_c[T]=V0pro_c/V0fee_c
}


#-------------------P11---------------------#
data<-data.frame(maturity<-seq(8),CDS_spread<-c(13.9,16.4,20.5,24,30.5,36,40.2,44.1)/10000)
colnames(data)<-c("maturity","CDS_spread")
rf=0.05
rec=0.4
L=1-rec
feetime <- seq(0.25,8,0.25)
maturity<-seq(1,8,1)
#------------------------
lamd0<-data$CDS_spread[1]/L
#for lamd1
t1<-seq(0.25,1,by=0.25)
V1fee<-sum(exp(-rf*t1)*exp(-lamd0*t1))
data$lamb<-calibrate.cds(rf, t=feetime, T=maturity, cdsrate=data$CDS_spread, RR=rec)





lambUse<-function(lamb0,lamb1,lamb2,lamb3,lamb4,lamb5,lamb6,lamb7,T){
  lambUse = ifelse(T<1,lamb0*T,ifelse(x<2,lamb0+lamb1*(T-1), 
                  ifelse(x<3,lamb0+lamb1+lamb2*(T-2),ifelse(x<4,lamb0+lamb1+lamb2+lamb3*(T-3),
                  ifelse(x<5,lamb0+lamb1+lamb2+lamb3+lamb4*(T-4),ifelse(x<6,lamb0+lamb1+lamb2+lamb3+lamb4+lamb5*(T-5),
                  ifelse(x<7,lamb0+lamb1+lamb2+lamb3+lamb4+lamb5+lamb6*(T-6),lamb0+lamb1+lamb2+lamb3+lamb4+lamb5+lamb6+lamb7*(T-7)                                                     
                  )))))))
  
}


feelamb<-c()
for (i in 1:32){feelamb<-lambUse(lamb0,lamb1,lamb2,lamb3,lamb4,lamb5,lamb6,lamb7,feetime(i))}
  
  
#---------------(c)-------------------#
#CIR model
CIRsprd<-function(kappa,theta,sigma,lambda0,steps,T){
  lambda_current = lambda0
  dt = 1/steps
  for (i in 1:T*steps){
      dwt = sqrt(dt)*rnorm(1,0,1)
      dlambda = kappa*(theta - lambda_current) + sigma* sqrt(lambda_current) * dwt
      lambda = lambda_current + dlambda
      lambda_current = lambda
  }
}










  
  


