### log ########################################################################################
## May 6
# Q.stat : changed Q1 and Q2
## May 7 ~ 10
# truncnorm: added to library(0.1), added truncation criteria(1.1) and replace rnorm with truncnorm(2.1)
# Wilson-Score CI: added to function(0.2.4), added to sampling loop statement(2.1)
# CI summary function: extract repeating part into functions(0.3)
# Coverage probability function: extract repeating part into functions(0.4)
# added exponential samples(2.1) => but something strange...
## May 21
# (1.1) typo - parameter of exp dist: c.x,y ... beta.x,y -> 1/beta.x.y
# (2.1) n.sim = 1000: for efficiency set n.sim=1000 instead of 10,000 until we finalize the code.
# (0.1) added package "knitr", "markdown"
## May 23
# (0.4.1) added labels on tables of coverage prob. function
# (0.4.2~3) added L/R NCP, ZWI function
## May 29
# (0.2.4~0.2.9) inserted CI.base at 0.2.4. renumbered from 0.2.4~0.2.8 to 0.2.5~0.2.9
# (0.2.10) added Clopper-Pearson at 0.2.10
# (0.2.4~0.3) revised overal functions
## June 10
# (0.2.11) added Bamber
# (0.2.12) added Wilson (with and without continuity-correction)
# (0.2.13) added Halperine Mee
# (0.2.14) added Double Beta
## June 16
# (0.2.9) RG. corrected z.a2 with z.a
# (0.2.15) added Double Gaussian


## 0. library  ########################################################################################
# 0.1 packages
# ROCR / dplyr
library(ROCR)
library(dplyr)
library(truncnorm)  # for 2.1 sampling
library(rootSolve)  # for Newton-Raphson Method (0.2.4)
library(sn)         # for Owen's T-function in Double Gaussian Method (0.2.15)

# 0.2 functions
# 0.2.1 rbivariate
rbivariate <- function(mu.x, sd.x, mu.y, sd.y, r=0.5, iter=100) {
  z1 <- rnorm(iter)
  z2 <- rnorm(iter)
  x <- sqrt(1-r^2)*sd.x*z1 + r*sd.x*z2 + mu.x
  y <- sd.y*z2 + mu.y
  return(list(x,y))
}
# 0.2.2 truncated exponential / qtexp: qunatile, rtexp: random sample of t.exp
# http://tolstoy.newcastle.edu.au/R/e10/help/10/06/8669.html
qtexp <- function(x, beta, trunc) { -log(1-x*(1-exp(-trunc/beta)))*beta }
rtexp <- function(n, beta, trunc) { qtexp(runif(n), beta, trunc) }

# 0.2.3 Q.stat : Q1 and Q2
# log: Error corrected on May 6: Q1 & Q2 switched
Q.stat <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker") {
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  Q2 <- Q1 <- 0
  for (i in 1:n.x) {
    Q1 = Q1 + ( sum(y > x[i]) + sum(y==x[i])/2 )^2
  }
  for (j in 1:n.y) {
    Q2 = Q2 + ( sum(y[j] > x) + sum(y[j]==x)/2 )^2
  }
  Q1 = Q1 /(n.x * n.y^2); Q2 = Q2 /(n.x^2 * n.y)
  return(data.frame(Q1 = Q1, Q2 = Q2))
}


# 0.2.4 Basic CI function
CI.base <- function(mu.hat, Var, alpha) {
  return(data.frame(lb = mu.hat - z.a2*sqrt(Var), ub = mu.hat + z.a2*sqrt(Var))) }

# 0.2.5 Wilson-Score CI
# V for variance, WS.equation for equation function, WS for CI bounds solution
V.WS <- function(AUC, n.x, n.y, Q1 = AUC/(2-AUC), Q2=2* AUC^2 /(AUC + 1)) {
  return((AUC * (1-AUC) + (n.y -1)*(Q1 - AUC^2) + (n.x-1)*(Q2 - AUC^2))/(n.x*n.y) )
}
WS.equation <- function(AUC, theta.hat, n.x, n.y, alpha) {
  AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.WS(AUC, n.x, n.y)) - theta.hat
}
WS <- function(theta.hat, n.x, n.y, alpha) {
  start.1 = max(0.01, theta.hat - .1); start.2 = min(0.99, theta.hat + .1)
  multiroot(WS.equation, c(start.1, start.2), theta.hat=theta.hat, n.x=n.x, n.y=n.y, alpha=alpha)$root 
}


# 0.2.6 Newcombe 1
V.NC1 <- function(AUC, n.x, n.y, Q1 = AUC/(2-AUC), Q2 = 2*AUC^2 /(AUC + 1), N = (n.x+n.y)/2) {
  return(max(AUC*(1-AUC)/(n.x*n.y)*((2*N-1)- (3*N-3)/((2-AUC)*(AUC+1))),0))
}


# 0.2.7 Newcombe 2
V.NC2 <- function(AUC, n.x, n.y, Q1 = AUC/(2-AUC), Q2=2* AUC^2 /(AUC + 1)) {
  return((AUC * (1-AUC) + (1/2*(n.x +n.y)-1)*(Q1 - AUC^2) + (1/2*(n.x +n.y)-1)*(Q2 - AUC^2))/(n.x*n.y) )
}
NC2.equation <- function(AUC, theta.hat, n.x, n.y, alpha) {
  AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.NC2(AUC, n.x, n.y)) - theta.hat   
}
NC2 <- function(theta.hat, n.x, n.y, alpha) {
  start.1 = max(0.01, theta.hat - .1); start.2 = min(0.99, theta.hat + .1)
  multiroot(NC2.equation, c(start.1, start.2), theta.hat=theta.hat, n.x=n.x, n.y=n.y, alpha=alpha)$root 
}


# 0.2.8 Cortes-Mohri
V.CM <- function(AUC, n.x, n.y){
  z1 = function (w, n.x, n.y, k) {1-0.5*(w/n.x+(k-w)/n.y)}
  z2 = function (w, n.x, n.y, k) 
    (n.y*w^2+n.x*(k-w)^2+n.y*(n.y+1)*w+n.x*(n.x+1)*(k-w)-2*w*(k-w)*(n.y+n.x+1))/(12*n.y^2*n.x^2)
  F.base <- function (w, n.x, n.y, k) {choose(n.y-k+2*w,w)*choose(n.x+k-2*w,k-w)}
  k=round((n.x+n.y)*(1-AUC))
  F1=F2=F3=F.den=0
  for (w in 0:k){
    F1 = F1 + F.base(w, n.x, n.y, k)*z1(w, n.x, n.y, k)^2
    F.den = F.den + F.base(w, n.x, n.y, k)
    F2 = F2 + F.base(w, n.x, n.y, k)*z1(w, n.x, n.y, k)
    F3 = F3 + F.base(w, n.x, n.y, k)*z2(w, n.x, n.y, k)
  }
  return(V = (F1/F.den - F2^2/F.den^2 + F3/F.den))
}


# 0.2.9 Reiser-Guttman    (corrected: z.a2 -> z.a)
RG <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker", alpha){
  z.a = qnorm(1-alpha)
  mu.x = mean(data[data[,disease]==0,marker])
  mu.y = mean(data[data[,disease]==1,marker])
  sig.x = sd(data[data[,disease]==0,marker])
  sig.y = sd(data[data[,disease]==1,marker])
  s2 = sig.x^2 + sig.y^2
  d = (mu.y-mu.x)/sqrt(s2)
  f = s2^2/((sig.x^4/(n.x-1))+(sig.y^4/(n.y-1)))
  M = s2/((sig.x^2/n.x)+(sig.y^2/n.y))
  delta1 = d - z.a*sqrt((1/M)+(d^2/(2*f)))
  delta2 = d + z.a*sqrt((1/M)+(d^2/(2*f)))
  return(data.frame(lb=pnorm(delta1), ub=pnorm(delta2)))
}


# 0.2.10 Clopper-Pearson
CP <- function (data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker", alpha){
  AUC.hat = (prediction(data[,marker], data[,disease]) %>% performance("auc"))@y.values[[1]]
  n = n.x + n.y; k = round(AUC.hat * n)
  f1 = qf(alpha/2,2*k,2*(n-k+1))
  f2 = qf(1-alpha/2,2*(k+1),2*(n-k))
  return(data.frame(lb=(k*f1)/(n-k+1+k*f1), ub=((k+1)*f2)/(n-k+(k+1)*f2)))
}


# 0.2.11 Bamber
V.Bm <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker",
                 AUC.hat = (prediction(data[,marker], data[,disease]) %>% performance("auc"))@y.values[[1]]){
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  b.yyx <- b.xxy <- p.xy <- 0
  for (i in 1:n.x) {b.yyx = b.yyx + sum(y < x[i])^2 + sum(x[i] < y)^2 -2*sum(y < x[i])*sum(x[i] < y) }
  b.yyx = b.yyx/(n.x*n.y^2)
  for (j in 1:n.y) {b.xxy = b.xxy + sum(x < y[j])^2 + sum(y[j] < x)^2 -2*sum(x < y[j])*sum(y[j] < x) }
  b.xxy = b.xxy/(n.y*n.x^2)
  for (i in 1:n.x) {  p.xy = p.xy + sum(y != x[i] ) }
  p.xy = p.xy /(n.x*n.y)
  return(1/(4*(n.y-1)*(n.x-1))*(p.xy + (n.y -1)*b.xxy + (n.x -1)*b.yyx - 4*(n.y + n.x -1)*(AUC.hat -0.5)^2))
}

# 0.2.12 Wilson
Wilson <- function(AUC.hat, n.x, n.y, alpha, cc=FALSE){
  n = n.x + n.y ; z.a2 = qnorm(1-alpha/2)
  t = z.a2^2 / n
  return(data.frame(lb = (2*n*AUC.hat + z.a2^2 - 1*cc - z.a2* sqrt(z.a2^2 - 2*cc - 1/n*cc + 4*AUC.hat*(n*(1-AUC.hat) +1*cc)))/(2*(n+z.a2^2)),
                    ub = (2*n*AUC.hat + z.a2^2 + 1*cc + z.a2* sqrt(z.a2^2 + 2*cc - 1/n*cc + 4*AUC.hat*(n*(1-AUC.hat) +1*cc)))/(2*(n+z.a2^2))))
}


# 0.2.13 Mee
# by definition.. something wrong with p1
p.stat2 <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker") {
  a <- Sys.time()
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  combn.y <- combn(n.y,2)
  combn.x <- combn(n.x,2)
  for (i in 1:(length(combn.x)/2)) {
    p1 = p1 + ((y > x[combn.x[1,i]]) + (y==x[combn.x[1,i]])/2 ) %*% ( (y > x[combn.x[2,i]]) + (y==x[combn.x[2,i]])/2)
  }
  for (i in 1:(length(combn.y)/2)) {
    p2 = p2 + ((x < y[combn.y[1,i]]) + (x==y[combn.y[1,i]])/2 ) %*% ( (x < y[combn.y[2,i]]) + (x==y[combn.y[2,i]])/2)
  }
  p1 = p1 *2 /(n.x * (n.x - 1) * n.y ); p2 = p2*2 /(n.x * n.y * (n.y - 1))
  print(Sys.time()-a)
  return(data.frame(p1 = p1, p2 = p2))
}
# by relation with U and Q statistics : much faster than p.stat2
p.stat <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker") {
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  Q = Q.stat(data,n.x,n.y,disease,marker)
  Q1 <- Q$Q1;  Q2 <- Q$Q2
  y <- matrix(rep(y,n.x),ncol=n.x)
  x <- t(matrix(rep(x,n.y),nrow=n.x))
  U.sq <- sum(((y>x)+(y==x)/2)^2)
  p2 <- (n.x*n.y^2*Q1 - U.sq)/(n.x*n.y*(n.y-1))
  p1 <- (n.y*n.x^2*Q2 - U.sq)/(n.y*n.x*(n.x-1))
  return(data.frame(p1=p1,p2=p2))
}

# by relation with U and Q statistics : much faster than p.stat2
Mee.stat <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker") {
  x = data[data[,disease]==0, marker]
  y = data[data[,disease]==1, marker]
  ranks = rank(c(x,y))
  rank.x <- ranks[1:n.x]
  rank.y <- ranks[-1:-n.x]
  AUC.hat <- (prediction(data[,marker], data[,disease]) %>% performance("auc"))@y.values[[1]]
  increment = ifelse(AUC.hat < 0.5, +0.5, -0.5)
  while (min(AUC.hat,1-AUC.hat)*sqrt(n.x*n.y) < 0.5){
    rank.y = rank.y + increment
    x <- data[data[,disease]==0, marker] <- rank.x
    y <- data[data[,disease]==1, marker] <- rank.y
    AUC.hat= (prediction(data[,marker], data[,disease]) %>% performance("auc"))@y.values[[1]]
  }
  p = p.stat(data, n.x, n.y)
  rho.hat = (p-AUC.hat^2)/(AUC.hat-AUC.hat^2)
  N.J.hat = n.x*n.y/(((n.x-1)*rho.hat$p1 +1)/(1-1/n.y) + ((n.y-1)*rho.hat$p2 + 1 )/(1-1/n.x) )
  return(data.frame(p1=p$p1, p2=p$p2, AUC.hat=AUC.hat, N.J.hat=N.J.hat))
}  
Mee.equation <- function(AUC, theta.hat, N.J.hat, alpha){
  AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(AUC*(1-AUC)/N.J.hat) - theta.hat
}
Mee <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker", alpha) {
  Mee.stat <- Mee.stat(data, n.x, n.y, disease, marker)
  theta.hat = Mee.stat$AUC.hat
  N.J.hat = Mee.stat$N.J.hat
  start.1 = max(0.01, theta.hat - .1); start.2 = min(0.99, theta.hat + .1)
  multiroot(Mee.equation, c(start.1, start.2), theta.hat=theta.hat, N.J.hat=N.J.hat, alpha=alpha)$root 
}


# 0.2.14 Double exponential

V.DB <- function(alp, n.x, n.y) {
  R1 = gamma(alp + 1)^2 / gamma(2*alp + 1)
  R2 = gamma(2*alp +1)*gamma(alp +1)/gamma(3*alp+1)
  theta = 1 - R1
  Q = 1 - 2*R1 + R2
  V = (theta*(1-theta) + (n.x+n.y-2)*(Q-theta^2))/n.x/n.y
}

DB.equation <- function(alp, theta.hat, n.x, n.y, alpha) {
  AUC = 1 - gamma(alp + 1)^2 / gamma(2*alp + 1)
  return(AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.DB(alp, n.x, n.y)) - theta.hat)
}
DB <- function(theta.hat, n.x, n.y, alpha) {
  start.1 = 1; start.2 = 3
  alp <- multiroot(DB.equation, c(start.1, start.2), theta.hat=theta.hat, n.x=n.x, n.y=n.y, alpha=alpha)$root
  return(1 - gamma(alp + 1)^2 / gamma(2*alp + 1))
}


# 0.2.15 Double Gaussian

V.DG <- function(delta, n.x, n.y) {
  theta = pnorm(delta/sqrt(2))
  Owen = T.Owen(delta/sqrt(2), 1/sqrt(3))
  return(V = ((n.x + n.y - 1)*theta*(1-theta) - 2*(n.x + n.y -2) * Owen ) /n.x / n.y)
}

DG.equation <- function(delta, theta.hat, n.x, n.y, alpha) {
  AUC = pnorm(delta/sqrt(2))
  return(AUC + c(+1,-1)*qnorm(1-alpha/2) * sqrt(V.DG(delta, n.x, n.y)) - theta.hat)
}

DG <- function(theta.hat, n.x, n.y, alpha) {
  start.1 = 1; start.2 = 3
  delta <- multiroot(DG.equation, c(start.1, start.2), theta.hat=theta.hat, n.x=n.x, n.y=n.y, alpha=alpha)$root
  return(pnorm(delta/sqrt(2)))
}



# 0.3 CI summary function  #################################################################################
CI.method=c("HMW", "HM", "WS", "NC1", "NC2", "CM", "RG", "CP", "Bm", "WS2", "WS3", "Mee", "DB", "DG")

AUC.CI <- function(data, n.x=sum(data[,disease]==0), n.y=sum(data[,disease]==1), disease="disease", marker="marker", alpha) {
  z.a2 = qnorm(1-alpha/2)
  AUC.hat = (prediction(data[,marker], data[,disease]) %>% performance("auc"))@y.values[[1]]  # package ROCR
  Q1 = Q.stat(data,n.x,n.y)$Q1 ; Q2 = Q.stat(data,n.x,n.y)$Q2                   # function Q.stat
  
  V.HMW = V.WS(AUC.hat, n.x, n.y, Q1, Q2)
  CI.HMW = CI.base(AUC.hat, V.HMW, alpha)
  
  V.HM = V.WS(AUC.hat, n.x, n.y)           # Q1, Q2 as default
  CI.HM = CI.base(AUC.hat, V.HM, alpha)
  
  CI.WS = data.frame(lb=rep(NA,3), ub=rep(NA,3)); for (z in 1:3) (CI.WS[z,] = WS(AUC.hat, n.x, n.y, alpha[z]))
  
  V.NC1 = V.NC1(AUC.hat, n.x, n.y)
  CI.NC1 = CI.base(AUC.hat, V.NC1, alpha)
  
  CI.NC2 = data.frame(lb=rep(NA,3), ub=rep(NA,3)); for (z in 1:3) (CI.NC2[z,] = NC2(AUC.hat, n.x, n.y, alpha[z]))
  
  V.CM = V.CM(AUC.hat, n.x, n.y)
  CI.CM = CI.base(AUC.hat, V.CM, alpha)
  
  CI.RG = RG(data, n.x, n.y, alpha=alpha) 
  
  CI.CP = CP(data, n.x, n.y, alpha=alpha)
  
  V.Bm = V.Bm(data, n.x, n.y, AUC.hat=AUC.hat)
  CI.Bm = CI.base(AUC.hat, V.Bm, alpha)
  
  CI.WS2 = Wilson(AUC.hat, n.x, n.y, alpha, F)
  CI.WS3 = Wilson(AUC.hat, n.x, n.y, alpha, T)
  
  CI.Mee = data.frame(lb=rep(NA,3), ub=rep(NA,3)); for (z in 1:3) (CI.Mee[z,] = Mee(data, n.x, n.y, alpha=alpha[z]))
  
  CI.DB = data.frame(lb=rep(NA,3), ub=rep(NA,3)); for (z in 1:3) (CI.DB[z,] = DB(AUC.hat, n.x, n.y, alpha=alpha[z]))
  
  CI.DG = data.frame(lb=rep(NA,3), ub=rep(NA,3)); for (z in 1:3) (CI.DG[z,] = DG(AUC.hat, n.x, n.y, alpha=alpha[z]))
  
  return(data.frame(AUC.hat=AUC.hat, n.x = n.x, n.y = n.y, Q1 = Q1 , Q2 = Q2 ,
                    HMW.lb = CI.HMW$lb,  HMW.ub = CI.HMW$ub,
                    HM.lb  = CI.HM$lb,   HM.ub = CI.HM$ub,
                    WS.lb  = CI.WS$lb,   WS.ub = CI.WS$ub,
                    NC1.lb = CI.NC1$lb,  NC1.ub = CI.NC1$ub,
                    NC2.lb = CI.NC2$lb,  NC2.ub = CI.NC2$ub,
                    CM.lb  = CI.CM$lb,   CM.ub = CI.CM$ub,
                    RG.lb  = CI.RG$lb,   RG.ub = CI.RG$ub,
                    CP.lb  = CI.CP$lb,   CP.ub = CI.CP$ub,
                    Bm.lb  = CI.Bm$lb,   Bm.ub = CI.Bm$ub,
                    WS2.lb = CI.WS2$lb,  WS2.ub = CI.WS2$ub,
                    WS3.lb = CI.WS3$lb,  WS3.ub = CI.WS3$ub,
                    Mee.lb = CI.Mee$lb,  Mee.ub = CI.Mee$ub,
                    DB.lb  = CI.DB$lb,   DB.ub = CI.DB$ub,
                    DG.lb  = CI.DG$lb,   DG.ub = CI.DG$ub
  ))
}

# 0.4 evaluation function  #################################################################################
# 0.4.1 CP(Coverage probability)
coverage <- function(data, dim = c(5,9,3), CI.method=CI.method) {
  # data: list of inferences(AUC hat, CI's lower/upper bounds)
  coverage <- array(,c(dim,length(CI.method)))
  dimnames(coverage) <- list(
    paste0("AUC = ",c(0.6, 0.7, 0.8, 0.9, 0.95)), 
    paste0(rep(c(30,50,200),each=3), "*", rep(c(0.3,0.5,0.7))),
    paste0("alpha = ",c("10%","5%","1%")),
    paste0("method = ",CI.method) )
  for (z in 1:length(CI.method)) {
    CI.lb = paste0(CI.method[z],".lb") ; CI.ub = paste0(CI.method[z],".ub")
    for (i in 1:dim[1]) {        # AUC (0.6 ~ 0.95)
      for (j in 1:dim[2]){       # n*pi (30,50,200) * (0.3, 0.5, 0.7)
        for (k in 1:dim[3]){     # alpha 10%, 5%, 1%
          coverage[i,j,k, z] <- mean(sapply(data[[i]][[j]],function(x) {
            x[k,CI.lb] <= x$AUC & x[k, CI.ub] >= x$AUC }), na.rm=T)
        }
      }
    }
  }
  cat(" row: AUC = 0.6 ~ 0.95 \n col: n*pi = (30,50,200)*(.3,.5,.7) \n 3rd: alpha = .1, .05, .01 \n 4th: CI methods = ", CI.method,  " \n\n ")
  return(coverage)
}

# 0.4.2 LRNCP(LNCP: left noncoverage prob, RNCP: right noncoverage prob)
LRNCP <- function(data, dim = c(5,9,3), CI.method=CI.method) {
  coverage <- array(,c(dim,length(CI.method)))
  dimnames(coverage) <- list(
    paste0("AUC = ",c(0.6, 0.7, 0.8, 0.9, 0.95)), 
    paste0(rep(c(30,50,200),each=3), "*", rep(c(0.3,0.5,0.7))),
    paste0("alpha = ",c("10%","5%","1%")),
    paste0("method = ",CI.method))
  for (z in 1:length(CI.method)) {
    CI.lb = paste0(CI.method[z],".lb") ; CI.ub = paste0(CI.method[z],".ub")
    for (i in 1:dim[1]) {        # AUC (0.6 ~ 0.95)
      for (j in 1:dim[2]){       # n*pi (30,50,200) * (0.3, 0.5, 0.7)
        for (k in 1:dim[3]){     # alpha 10%, 5%, 1%
          LNCP <- mean(sapply(data[[i]][[j]],function(x) {x[k,CI.lb] > x$AUC }), na.rm=T) %>% round(4) %>% format(nsmall=4)
          RNCP <- mean(sapply(data[[i]][[j]],function(x) {x[k,CI.ub] < x$AUC }), na.rm=T) %>% round(4) %>% format(nsmall=4)
          coverage[i,j,k, z] <- paste0(sub("^(-?)0.", "\\1.",LNCP),"+",sub("^(-?)0.", "\\1.",RNCP))
        }
      }
    }
  }
  cat(" row: AUC = 0.6 ~ 0.95 \n col: n*pi = (30,50,200)*(.3,.5,.7) \n 3rd: alpha = .1, .05, .01 \n 4th: CI methods = ", CI.method, " \n\n ")
  return(coverage)
}


# 0.4.3 ZWI rate(rate of occurence of zero width intervals)
ZWI <- function(data, dim = c(5,9,3), CI.method=CI.method) {
  zwi <- array(,c(dim,length(CI.method)))
  dimnames(zwi) <- list(
    paste0("AUC = ",c(0.6, 0.7, 0.8, 0.9, 0.95)), 
    paste0(rep(c(30,50,200),each=3), "*", rep(c(0.3,0.5,0.7))),
    paste0("alpha = ",c("10%","5%","1%")),
    paste0("method = ",CI.method))
  for (z in 1:length(CI.method)) {
    CI.lb = paste0(CI.method[z],".lb") ; CI.ub = paste0(CI.method[z],".ub")
    for (i in 1:dim[1]) {        # AUC (0.6 ~ 0.95)
      for (j in 1:dim[2]){       # n*pi (30,50,200) * (0.3, 0.5, 0.7)
        for (k in 1:dim[3]){     # alpha 10%, 5%, 1%
          zwi[i,j,k, z] <- mean(sapply(data[[i]][[j]],function(x) {x[k,CI.lb] == x[k,CI.ub] }), na.rm=T) %>% round(5) %>% format(nsmall=5)
        }
      }
    }
  }
  cat(" row: AUC = 0.6 ~ 0.95 \n col: n*pi = (30,50,200)*(.3,.5,.7) \n 3rd: alpha = .1, .05, .01 \n 4th: CI methods = ", CI.method, " \n\n ")
  return(zwi)
}






## 1. population data  ######################################################################################################
## 1.1 ~1.2 parameters
## 1.1 parms1 (theta, mu, sig, beta)        : dim(parm) = 5AUC + z.a2*sqrt(V.HM),

# theta = AUC
parm1 <- data.frame(theta = c(0.6, 0.7, 0.8, 0.9, 0.95))
parm1 <- within(parm1,
{ # Under normal distribution
  mu.x <- 1; sig.x <- 1; sig.y <- 2
  # Under exponential distribution
  beta.x <- 1
  # Residual params under thetas
  mu.y <- mu.x + qnorm(theta)*sqrt(sig.x^2 + sig.y^2)
  beta.y <- theta * beta.x / (1-theta)
  # bounds for trunction
  a.x <- mu.x - 3* sig.x; b.x <- mu.x + 3*sig.x   # bounds for truncated noraml
  a.y <- mu.y - 3* sig.y; b.y <- mu.y + 3*sig.y
  c.x <- qexp(0.99, 1/beta.x); c.y <- qexp(0.99, 1/beta.y) # bounds for truncated exp  # correction: beta. -> 1/beta. (scale parameter)
})

## 1.2 parms2: n, pi
# n= sample size,  pi = prevalence rate     : dim(parm2) = 3*3=9
# n.x and n.y (= n * pi) are assumed to be fixed not random.
parm2 <- data.frame(n = rep(c(30, 50, 200), each=3), pi = rep(c(0.30, 0.50, 0.70),3))


## 1.3 alpha levels
alpha <- c(0.1, 0.05, 0.01)
z.a2<- qnorm(1-alpha/2)     # z_alpha/2






## 2. Sample data  ##################################################################################################
## 2.1 generating sample data with bivariate-normal distribution

sample.normal <- sample.normal.stat <- list()   # outer shell   (data & statistics) normal
temp.2 <- temp.2.stat <- list()                 # middle shells
temp.1 <- temp.1.stat <- list()                 # inner shells
sample.exp <- sample.exp.stat <- list()         # outer shell   (data & statistics) exponential
temp.2.exp <- temp.2.exp.stat <- list()         # middle shells
temp.1.exp <- temp.1.exp.stat <- list()         # inner shells
{
  set.seed(10)
  for (i in 1:5) { # i: AUC
    for (j in 1:9) { # j: n sample size, pi prevalence rate
      mu.x <- parm1$mu.x[i]; mu.y <- parm1$mu.y[i]; sig.x <- parm1$sig.x[i]; sig.y <- parm1$sig.y[i];
      beta.x <- parm1$beta.x[i]; beta.y <- parm1$beta.y[i]; theta <- parm1$theta[i]
      a.x <- parm1$a.x[i]; b.x <- parm1$b.x[i]; a.y <- parm1$a.y[i]; b.y <- parm1$b.y[i]
      c.x <- parm1$c.x[i]; c.y <- parm1$c.y[i]
      n <- parm2$n[j] ; pi <- parm2$pi[j]
      n.y <- n*pi; n.x <- n - n.y
      for (k in 1:1000){ # k: num of simulations(samples)
        # normal
        temp.1[[k]] <- temp <- data.frame(disease = c(rep(0,n.x), rep(1, n.y)),
                                          marker  = c(rtruncnorm(n-n*pi,a.x, b.x, mu.x, sig.x), rtruncnorm(n*pi,a.y, b.y, mu.y, sig.y)))
        temp.1.stat[[k]] <- AUC.CI(temp, n.x, n.y, alpha=alpha)
        temp.1.stat[[k]]$AUC <- theta
        # exponential
        temp.1.exp[[k]] <- temp <- data.frame(disease = c(rep(0,n.x), rep(1, n.y)),
                                              marker  = c(rtexp(n-n*pi, beta.x, c.x), rtexp(n*pi,beta.y, c.y)))
        temp.1.exp.stat[[k]] <- AUC.CI(temp, n.x, n.y, alpha=alpha)
        temp.1.exp.stat[[k]]$AUC <- theta
      }
      temp.2[[j]] <- temp.1
      temp.2.stat[[j]] <- temp.1.stat
      temp.2.exp[[j]] <- temp.1.exp
      temp.2.exp.stat[[j]] <- temp.1.exp.stat
    }
    sample.normal[[i]] <- temp.2
    sample.normal.stat[[i]] <- temp.2.stat
    sample.exp[[i]] <- temp.2.exp
    sample.exp.stat[[i]] <- temp.2.exp.stat
  }
}

rm(temp.2, temp.2.stat, temp.1, temp.1.stat, temp.2.exp, temp.2.exp.stat, temp.1.exp, temp.1.exp.stat) #clearing temporary data

## 3. Evaluation  ########################################################################################
#3.1 Coverage probability
print(coverage.normal <- coverage(sample.normal.stat, CI.method=CI.method))
print(coverage.exp <- coverage(sample.exp.stat, CI.method=CI.method))

#3.2 Left/Right Non Coverage Probability
print(LRNCP.normal <- LRNCP(sample.normal.stat, CI.method=CI.method))
print(LRNCP.exp <- LRNCP(sample.exp.stat, CI.method=CI.method))

#3.3 Rate of occurrence of Zero Width Intervals
print(ZWI.normal <- ZWI(sample.normal.stat, CI.method=CI.method))
print(ZWI.exp <- ZWI(sample.exp.stat, CI.method=CI.method))

