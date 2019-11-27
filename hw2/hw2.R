#question 1-------------------------------
r1 = 102/99-1

p2 = 98
f2 = function(x){
  2.5/(1+r1)+102.5/(1+x)^2 - p2
}

log((98-2.5*exp(-0.02985))/102.5)*(-0.5)
log((97.5-3*exp(-0.02985)-3*exp(-2*0.03498))/103)*(-1/3)
4*exp(-0.02985)+4*exp(-2*0.03498)+104*exp(-3*0.03841)

p = 97.5
f3 = function(x){
  3/(1+r1)+3/(1+r2)^2+103/(1+x)^3 - p
}

library(rootSolve)
r2 = multiroot(f = f2, start = 0)$root
r3 = multiroot(f=f3,start = 0)$root

library(ggplot2)
plot(x = c(1,2,3),y = c(r1,r2,r3), type = "o")
qplot(x = c(1,2,3),y= c(r1,r2,r3), geom = "line")

p4 = 4/(1+r1)+4/(1+r2)^2 + 104/(1+r3)^3

#question 2-------------------------------
liquid_r1 = 101/99-1

log(99/101)*(-1)
log((98-2*exp(0.02))/102)*(-0.5)
log((98.5-3*exp(0.02)-3*exp(-2*0.031))/103)*(-1/3)

liquid_p2 = 98
f2 = function(x){
  2/(1+liquid_r1)+102/(1+x)^2 - liquid_p2
}
liquid_r2 = multiroot(f=f2,start = 0)$root

liquid_p3 = 98.5
f3 = function(x){
  3/(1+liquid_r1)+3/(1+liquid_r2)^2+103/(1+x)^3 - liquid_p3
}
liquid_r3 = multiroot(f=f3,start = 0)$root

qplot(x = c(1,2,3),y= c(liquid_r1,liquid_r2,liquid_r3), geom = "line")

il_r1 = log(98/101)*(-1)

il_p2 = 97
f2 = function(x){
  2/exp(il_r1)+102/exp(2*x) - il_p2
}
il_r2 = multiroot(f=f2,start = 0)$root

il_p3 = 96.5
f3 = function(x){
  3/exp(il_r1)+3/exp(2*il_r2)+103/exp(3*x) - il_p3
}
il_r3 = multiroot(f=f3,start = 0)$root

4*exp(-0.02)+4*exp(-2*0.031)+104*exp(-3*0.03541)
4*exp(-il_r1)+4*exp(-2*il_r2)+104*exp(-3*il_r3)

qplot(x = c(1,2,3),y= c(il_r1,il_r2,il_r3), geom = "line")

f = function(c,r1,r2,r3){
  c/(1+r1)+c/(1+r2)^2+(100+c)/(1+r3)^3
}
liquid_p4 = f(4,liquid_r1,liquid_r2,liquid_r3)
il_p4 = f(4,il_r1,il_r2,il_r3)

#question 3-------------------------------
d1 = (log(140/100)+0.5*(0.2^2)*9)/(0.2*sqrt(9))
d2 = d1 - 0.2*sqrt(9)
N_nd1 = pnorm(-d1,mean = 0, sd = 1)
N_d2 = pnorm(d2,mean = 0, sd = 1)
Debt = 140*N_nd1+100*N_d2
#(b)
spread = (log(100)-log(Debt))/9
#(c)
N_d1 = pnorm(d1,mean = 0,sd = 1)
Equity = 140*N_d1-100*N_d2
#(d)
sigma_equity = N_d1*0.2*140/Equity

#question 4-------------------------------
r = log(95/100)*(-1)
R = (100+r)
u = 120/118
d = 60/118
pd = (R-d)/(u-d)
pu = 1-pd
