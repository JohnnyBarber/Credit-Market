#####################################################################################
#                                      part a                                       #
#####################################################################################
library(data.table)

compute_spread = function(time,intensity,R,r){
  l = 1 - R
  fixed = seq(0.25,1,0.25)
  t = seq(0.25,time,0.25)
  
  if(time < 0){
    stop("Error: number of years has to be positive!")
    
  }else if(length(intensity) > 8){
    stop("Error: exceed maximum steps of default!")
    
  }else if(time >= 0 && time <= 1){
    S_01 = l*intensity[1]
    return(S_01)
    
  }else if(time > 1 && time <= 2){
    S12 = function(x){
      Vf = sum(exp(-r*fixed)*exp(-x[1]*fixed)) + sum(exp(-r*(t+1))*exp(-x[2]-x*t))
      Vp = l*x[1]*sum(exp(-r*fixed)*exp(-x[1]*fixed)) +
        l*x[2]*sum(exp(-r*(t+1))*exp(-x[1]-x[2]*t))
      return(Vp/Vf)
    }
    S_12 = S12(intensity)
    return(S_12)
    
  }else if(time > 2 && time <= 3){
    S23 = function(x){
      Vf = sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed)) +
        sum(exp(-r*(t+2))*exp(-x[1]-x[2]-x[3]*t))
      Vp = l*x[1]*sum(exp(-r*fixed)*exp(-x[1]*fixed)) +
        l*x[2]*sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed))+
        l*x[3]*sum(exp(-r*(t+2))*exp(-x[1]-x[2]-x[3]*t))
      return(Vp/Vf)
    }
    S_23 = S23(intensity)
    return(S_23)
    
  }else if(time > 3 && time <= 4){
    S34 = function(x){
      Vf = sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed)) +
        sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        sum(exp(-r*(t+3))*exp(-x[1]-x[2]-x[3]-x[4]*t))
      Vp = l*x[1]*sum(exp(-r*fixed)*exp(-x[1]*fixed)) +
        l*x[2]*sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed))+
        l*x[3]*sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        l*x[4]*sum(exp(-r*(t+3))*exp(-x[1]-x[2]-x[3]-x[4]*t))
      return(Vp/Vf)
    }
    S_34 = S34(intensity)
    return(S_34)
    
  }else if(time > 4 && time <= 5){
    S45 = function(x){
      Vf = sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        sum(exp(-r*(fixed+1))*exp(-x[1]-x[1]*fixed)) +
        sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        sum(exp(-r*(t+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*t))
      Vp = l*x[1]*sum(exp(-r*fixed)*exp(-x[1]*fixed)) +
        l*x[2]*sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed))+
        l*x[3]*sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        l*x[4]*sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        l*x[5]*sum(exp(-r*(t+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*t))
      return(Vp/Vf)
    }
    S_45 = S45(intensity)
    return(S_45)
    
  }else if(time > 5 && time <= 6){
    S56 = function(x){
      Vf = sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed)) +
        sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        sum(exp(-r*(fixed+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*fixed))+
        sum(exp(-r*(t+5))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]*t))
      Vp = l*x[1]*sum(exp(-r*fixed)*exp(-x[1]*fixed)) +
        l*x[2]*sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed))+
        l*x[3]*sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        l*x[4]*sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        l*x[5]*sum(exp(-r*(fixed+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*fixed))+
        l*x[6]*sum(exp(-r*(t+5))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]*t))
      return(Vp/Vf)
    }
    S_56 = S56(intensity)
    return(S_56)
    
  }else if(time > 6 && time <= 7){
    S67 = function(x){
      Vf = sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed)) +
        sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        sum(exp(-r*(fixed+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*fixed))+
        sum(exp(-r*(fixed+5))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]*fixed))+
        sum(exp(-r*(t+6))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]-x[7]*t))
      Vp = l*x[1]*sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        l*x[2]*sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed)) +
        l*x[3]*sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        l*x[4]*sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        l*x[5]*sum(exp(-r*(fixed+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*fixed))+
        l*x[6]*sum(exp(-r*(fixed+5))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]*fixed))+
        l*x[7]*sum(exp(-r*(t+6))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]-x[7]*t))
      return(Vp/Vf)
    }
    S_67 = S67(intensity)
    return(S_67)
    
  }else if(time > 7 && time <= 8){
    S78 = function(x){
      Vf = sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed)) +
        sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        sum(exp(-r*(fixed+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*fixed))+
        sum(exp(-r*(fixed+5))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]*fixed))+
        sum(exp(-r*(fixed+6))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]-x[7]*fixed))+
        sum(exp(-r*(t+7))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]-x[7]-x[8]*t))
      Vp = l*x[1]*sum(exp(-r*fixed)*exp(-x[1]*fixed)) + 
        l*x[2]*sum(exp(-r*(fixed+1))*exp(-x[1]-x[2]*fixed)) +
        l*x[3]*sum(exp(-r*(fixed+2))*exp(-x[1]-x[2]-x[3]*fixed))+
        l*x[4]*sum(exp(-r*(fixed+3))*exp(-x[1]-x[2]-x[3]-x[4]*fixed))+
        l*x[5]*sum(exp(-r*(fixed+4))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]*fixed))+
        l*x[6]*sum(exp(-r*(fixed+5))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]*fixed))+
        l*x[7]*sum(exp(-r*(fixed+6))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]-x[7]*fixed))+
        l*x[8]*sum(exp(-r*(t+7))*exp(-x[1]-x[2]-x[3]-x[4]-x[5]-x[6]-x[7]-x[8]*t))
      return(Vp/Vf)
    }
    S_78 = S78(intensity)
    return(S_78)
    
  }else{
    stop("Error: exceed maximum number of years!")
  }
}

compute_spread(time = 8,c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08), R = 0.4, r = 0.05)

#####################################################################################
#                                      part b                                       #
#####################################################################################
maturity = seq(8)
CDS_spread = c(13.9,16.4,20.5,24,30.5,36,40.2,44.1)/10000
time = seq(0.25,8,0.25)
r = 0.05
Recovery = 0.4

#---------------------------------Auto Calculation----------------------------------#
library(CreditRisk)

lamb = calibrate.cds(r=r,t = time, T = maturity, cdsrate = CDS_spread, RR = Recovery)

period = c(0,1,2,3,4,4,5,6,7,8,8,9,10,11,12,12,13,14,15,16,16,17,18,
           19,20,20,21,22,23,24,24,25,26,27,28,28,29,30,31,32)
lambda = rep(lamb$min_int,each = 5)
plot(period, lambda, cex = 0.1, 
     xlab = "time period", ylab = "default intensity",
     main = "Intensity Process")
grid()
lines(period[1:5],lambda[1:5],col = "blue",lwd=2)
lines(period[6:10],lambda[6:10],col = "blue",lwd=2)
lines(period[11:15],lambda[11:15],col = "blue",lwd=2)
lines(period[16:20],lambda[16:20],col = "blue",lwd=2)
lines(period[21:25],lambda[21:25],col = "blue",lwd=2)
lines(period[26:30],lambda[26:30],col = "blue",lwd=2)
lines(period[31:35],lambda[31:35],col = "blue",lwd=2)
lines(period[36:40],lambda[36:40],col = "blue",lwd=2)

#---------------------------------Manual Calculation----------------------------------#
maturity = seq(8)
CDS_spread = c(13.9,16.4,20.5,24,30.5,36,40.2,44.1)/10000
time = seq(0.25,8,0.25)
r = 0.05
Recovery = 0.4
l = 1 - Recovery

compute_intensity = function(l,CDS_spread,r){
  #compute lambda 0
  lamb0 = CDS_spread[1]/l
  
  #compute lambda 1
  S12 = function(x){
    t = c(0.25,0.5,0.75,1)
    Vf = sum(exp(-r*t)*exp(-lamb0*t)) + sum(exp(-r*(t+1))*exp(-lamb0-x*t))
    Vp = l*lamb0*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0*s)),lowerLimit = 0, upperLimit = 1)$integral +
      l*x*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-x*(s-1))),lowerLimit = 1, upperLimit = 2)$integral
    Vp/(0.25*Vf)-CDS_spread[2]
  }
  
  lamb1 = uniroot(S12,c(0,1), extendInt = "yes")$root
  
  #compute lambda 2
  S23 = function(x){
    t = c(0.25,0.5,0.75,1)
    Vf = sum(exp(-r*t)*exp(-lamb0*t)) + 
      sum(exp(-r*(t+1))*exp(-lamb0-lamb1*t)) +
      sum(exp(-r*(t+2))*exp(-lamb0-lamb1-x*t))
    Vp = l*lamb0*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0*s)),lowerLimit = 0, upperLimit = 1)$integral +
      l*lamb1*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1*(s-1))),lowerLimit = 1, upperLimit = 2)$integral+
      l*x*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-x*(s-2))),lowerLimit = 2, upperLimit = 3)$integral
    Vp/(0.25*Vf) - CDS_spread[3]
  }
  
  lamb2 = uniroot(S23,c(0,1), extendInt = "yes")$root
  
  #compute lambda 3
  S34 = function(x){
    t = c(0.25,0.5,0.75,1)
    Vf = sum(exp(-r*t)*exp(-lamb0*t)) + 
      sum(exp(-r*(t+1))*exp(-lamb0-lamb1*t)) +
      sum(exp(-r*(t+2))*exp(-lamb0-lamb1-lamb2*t))+
      sum(exp(-r*(t+3))*exp(-lamb0-lamb1-lamb2-x*t))
    Vp = l*lamb0*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0*s)),lowerLimit = 0, upperLimit = 1)$integral +
      l*lamb1*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1*(s-1))),lowerLimit = 1, upperLimit = 2)$integral +
      l*lamb2*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2*(s-2))),lowerLimit = 2, upperLimit = 3)$integral +
      l*x*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-x*(s-3))),lowerLimit = 3, upperLimit = 4)$integral
    Vp/(0.25*Vf) - CDS_spread[4]
  }
  
  lamb3 = uniroot(S34,c(0,1), extendInt = "yes")$root
  
  #compute lambda 4
  S45 = function(x){
    t = c(0.25,0.5,0.75,1)
    Vf = sum(exp(-r*t)*exp(-lamb0*t)) + 
      sum(exp(-r*(t+1))*exp(-lamb0-lamb1*t)) +
      sum(exp(-r*(t+2))*exp(-lamb0-lamb1-lamb2*t))+
      sum(exp(-r*(t+3))*exp(-lamb0-lamb1-lamb2-lamb3*t))+
      sum(exp(-r*(t+4))*exp(-lamb0-lamb1-lamb2-lamb3-x*t))
    Vp = l*lamb0*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0*s)),lowerLimit = 0, upperLimit = 1)$integral +
      l*lamb1*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1*(s-1))),lowerLimit = 1, upperLimit = 2)$integral+
      l*lamb2*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2*(s-2))),lowerLimit = 2, upperLimit = 3)$integral+
      l*lamb3*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3*(s-3))),lowerLimit = 3, upperLimit = 4)$integral+
      l*x*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-x*(s-4))),lowerLimit = 4, upperLimit = 5)$integral
    Vp/(0.25*Vf) - CDS_spread[5]
  }
  
  lamb4 = uniroot(S45,c(0,1), extendInt = "yes")$root
  
  #compute lambda 5
  S56 = function(x){
    t = c(0.25,0.5,0.75,1)
    Vf = sum(exp(-r*t)*exp(-lamb0*t)) + 
      sum(exp(-r*(t+1))*exp(-lamb0-lamb1*t)) +
      sum(exp(-r*(t+2))*exp(-lamb0-lamb1-lamb2*t))+
      sum(exp(-r*(t+3))*exp(-lamb0-lamb1-lamb2-lamb3*t))+
      sum(exp(-r*(t+4))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4*t))+
      sum(exp(-r*(t+5))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-x*t))
    Vp = l*lamb0*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0*s)),lowerLimit = 0, upperLimit = 1)$integral +
      l*lamb1*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1*(s-1))),lowerLimit = 1, upperLimit = 2)$integral+
      l*lamb2*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2*(s-2))),lowerLimit = 2, upperLimit = 3)$integral+
      l*lamb3*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3*(s-3))),lowerLimit = 3, upperLimit = 4)$integral+
      l*lamb4*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4*(s-4))),lowerLimit = 4, upperLimit = 5)$integral+
      l*x*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-x*(s-5))),lowerLimit = 5, upperLimit = 6)$integral
    Vp/(0.25*Vf) - CDS_spread[6]
  }
  
  lamb5 = uniroot(S56,c(0,1), extendInt = "yes")$root
  
  #compute lambda 6
  S67 = function(x){
    t = c(0.25,0.5,0.75,1)
    Vf = sum(exp(-r*t)*exp(-lamb0*t)) + 
      sum(exp(-r*(t+1))*exp(-lamb0-lamb1*t)) +
      sum(exp(-r*(t+2))*exp(-lamb0-lamb1-lamb2*t))+
      sum(exp(-r*(t+3))*exp(-lamb0-lamb1-lamb2-lamb3*t))+
      sum(exp(-r*(t+4))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4*t))+
      sum(exp(-r*(t+5))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5*t))+
      sum(exp(-r*(t+6))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5-x*t))
    Vp = l*lamb0*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0*s)),lowerLimit = 0, upperLimit = 1)$integral +
      l*lamb1*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1*(s-1))),lowerLimit = 1, upperLimit = 2)$integral+
      l*lamb2*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2*(s-2))),lowerLimit = 2, upperLimit = 3)$integral+
      l*lamb3*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3*(s-3))),lowerLimit = 3, upperLimit = 4)$integral+
      l*lamb4*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4*(s-4))),lowerLimit = 4, upperLimit = 5)$integral+
      l*lamb5*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5*(s-5))),lowerLimit = 5, upperLimit = 6)$integral+
      l*x*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5-x*(s-6))),lowerLimit = 6, upperLimit = 7)$integral
    Vp/(0.25*Vf) - CDS_spread[7]
  }
  
  lamb6 = uniroot(S67,c(0,1), extendInt = "yes")$root
  
  #compute lambda 7
  S78 = function(x){
    t = c(0.25,0.5,0.75,1)
    Vf = sum(exp(-r*t)*exp(-lamb0*t)) + 
      sum(exp(-r*(t+1))*exp(-lamb0-lamb1*t)) +
      sum(exp(-r*(t+2))*exp(-lamb0-lamb1-lamb2*t))+
      sum(exp(-r*(t+3))*exp(-lamb0-lamb1-lamb2-lamb3*t))+
      sum(exp(-r*(t+4))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4*t))+
      sum(exp(-r*(t+5))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5*t))+
      sum(exp(-r*(t+6))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5-lamb6*t))+
      sum(exp(-r*(t+7))*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5-lamb6-x*t))
    Vp = l*lamb0*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0*s)),lowerLimit = 0, upperLimit = 1)$integral +
      l*lamb1*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1*(s-1))),lowerLimit = 1, upperLimit = 2)$integral+
      l*lamb2*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2*(s-2))),lowerLimit = 2, upperLimit = 3)$integral+
      l*lamb3*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3*(s-3))),lowerLimit = 3, upperLimit = 4)$integral+
      l*lamb4*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4*(s-4))),lowerLimit = 4, upperLimit = 5)$integral+
      l*lamb5*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5*(s-5))),lowerLimit = 5, upperLimit = 6)$integral+
      l*lamb6*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5-lamb6*(s-6))),lowerLimit = 6, upperLimit = 7)$integral+
      l*x*adaptIntegrate(function(s)(exp(-r*s)*exp(-lamb0-lamb1-lamb2-lamb3-lamb4-lamb5-lamb6-x*(s-7))),lowerLimit = 7, upperLimit = 8)$integral
    Vp/(0.25*Vf) - CDS_spread[8]
  }
  
  lamb7 = uniroot(S78,c(0,1), extendInt = "yes")$root
  
  lambda_t = c(lamb0,lamb1,lamb2,lamb3,lamb4,lamb5,lamb6,lamb7)
  
  return(lambda_t)
}

plot_intensity = function(lambda){
  plt_lambda = rep(lambda, each = 5)
  period = c(0,1,2,3,4,4,5,6,7,8,8,9,10,11,12,12,13,14,15,16,16,17,18,
             19,20,20,21,22,23,24,24,25,26,27,28,28,29,30,31,32)
  plot(period, plt_lambda, cex = 0.1, 
       xlab = "time period", ylab = "default intensity",
       main = "Intensity Process")
  grid()
  lines(period[1:5],plt_lambda[1:5],col = "blue",lwd=2)
  lines(period[6:10],plt_lambda[6:10],col = "blue",lwd=2)
  lines(period[11:15],plt_lambda[11:15],col = "blue",lwd=2)
  lines(period[16:20],plt_lambda[16:20],col = "blue",lwd=2)
  lines(period[21:25],plt_lambda[21:25],col = "blue",lwd=2)
  lines(period[26:30],plt_lambda[26:30],col = "blue",lwd=2)
  lines(period[31:35],plt_lambda[31:35],col = "blue",lwd=2)
  lines(period[36:40],plt_lambda[36:40],col = "blue",lwd=2)
}

# lamb = c(CDS_spread[1]/l)
# for(i in 1:8){
#   Vf = 0
#   Vp = 0
#   compute = function(x){
#     for(j in 1:i){
#       Vf = Vf + exp(-lamb[j])*(exp(-r*j)/r-exp(-r*(j+1))/r)
#       Vp = Vp + l*lamb[j]*exp(-lamb[j])*(exp(-r*j)/r-exp(-r*(j+1))/r)
#     }
#     Vf = Vf + exp(-x)*(exp(-r*i)/r-exp(-r*(i+1))/r)
#     Vp = Vp + l*x*exp(-x)*(exp(-r*i)/r-exp(-r*(i+1))/r)
#     return(Vp/Vf-CDS_spread[i])
#   }
#   lamb[i] = uniroot(compute,c(0,1))$root
# }


lambda = compute_intensity(l,CDS_spread,r)
plot_intensity(lambda)

test_spread = seq(8)
for(i in 1:8){
  test_spread[i] = compute_spread(time = i,lambda, R = 0.4, r = 0.05)
}

#compare with the formula in part a
compare = data.table(CDS_spread,test_spread)
compare[, difference := abs(CDS_spread-test_spread)]
plot(seq(length(CDS_spread)),test_spread, type = "l",
     xlab = "time of maturity", ylab = "spread",
     main = "CDS Spread")
grid()
lines(seq(length(CDS_spread)),CDS_spread, col = "green")
lines(seq(length(CDS_spread)),test_spread, col = "red")
legend("topleft",legend = c("Actual","Backtested"),
       lty = 1, col = c("green","red"))

#####################################################################################
#                                      part c                                       #
#####################################################################################

