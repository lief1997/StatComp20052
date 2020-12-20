## -----------------------------------------------------------------------------
x <- rnorm(10)
y <- rnorm(10)
plot(x,y)

## -----------------------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
summary(lm.D9)$coef

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
x <- 2/sqrt(1-u)#F(x)=1-(2/x)^2
hist(x,prob=TRUE,breaks=100,xlim=c(0,20),main=expression(f(x)==8/x^3))#cut the hist in the range from 0 to 40
y <- seq(0,20,0.1)
lines(y, 8/y^3)



## -----------------------------------------------------------------------------
gen <- function(u1,u2,u3)
{
  u <- vector(length = length(u1))
  for(i in 1:length(u1)){
    if(abs(u3[i]) >= abs(u2[i]) & abs(u3[i]) >= abs(u1[i]))
     u[i] = u2[i] else
     u[i] = u3[i]
  }
  return(u)
}


## -----------------------------------------------------------------------------
n <- 1000
u1 <- runif(n, min = -1, max = 1)
u2 <- runif(n, min = -1, max = 1)
u3 <- runif(n, min = -1, max = 1)
u <- gen(u1,u2,u3)
hist(u, prob = TRUE)
y <- seq(-1, 1, .01)
lines(y, 3/4*(1-y^2))

## -----------------------------------------------------------------------------
m <- 1e4
x <- runif(m, min=0, max= pi/3)
theta.hat <- mean(sin(x)) * pi/3
print(c(theta.hat,cos(0) - cos(pi/3)))

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 1e4/2; x <- runif(m, min=0, max=1)
MC1 <- mean(exp(x)+exp(1-x))/2
print(MC1)

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 1e4; x <- runif(m, min=0, max=1)
MC2 <- mean(exp(x))
print(MC2)

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 1e4; x1 <- runif(m/2, min=0, max=1)
x2 <- runif(m, min=0, max=1)
print(1-var((exp(x1)+exp(1-x1))/2)/var(exp(x2)))

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- se <- numeric(2) 
g <- function(x) {
  x^2*exp(-x^2/2)/sqrt(2*pi) * (x>1) }

set.seed(1234)

x <- rexp(m, 1) 
fg <- g(x) / exp(-x) 
theta.hat[1] <- mean(fg) 
se[1] <- sd(fg)

x <- rnorm(m)
fg <- g(x) / (exp(-x^2/2)/sqrt(2*pi))
theta.hat[2] <- mean(fg) 
se[2] <- sd(fg)

rbind(theta.hat, se)

## -----------------------------------------------------------------------------
m <- 10000
N <- 5
T <- matrix(0,N,2)
for(i in 1:N)
{
  u <- runif(m)
  z <- exp((1-i)/5)-exp(-i/5) 
  x <- -log(exp(-0.2*i+0.2)-(u*z)) 
  g <- function(y){z/(1+y^2)*(y>((i-1)/5))*(y<(i/5))}
  d <-g(x)
  T[i,1] <- mean(d)
  T[i,2] <- var(d)
}
est1 <- sum(T[,1]);est1
est2 <- 5*sum(T[,2]);est2


## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
set.seed(1234)
UCL <- replicate(1000, expr = {x <- rnorm(n, mean = 0, sd = 1)
mean(x) + var(x)*qt(alpha/2, df = n-1)/sqrt(n) })
DCL <- replicate(1000, expr = {x <- rnorm(n, mean = 0, sd = 1)
mean(x) - var(x)*qt(alpha/2, df = n-1)/sqrt(n) })
mean(UCL < 0 & 0 < DCL)


## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
set.seed(1234)
UCL <- replicate(1000, expr = {x <- rchisq(n, df=2)
mean(x) + var(x)*qt(alpha/2, df = n-1)/sqrt(n) })
DCL <- replicate(1000, expr = {x <- rchisq(n, df=2)
mean(x) - var(x)*qt(alpha/2, df = n-1)/sqrt(n) })
mean(UCL < 0 & 0 < DCL)

## -----------------------------------------------------------------------------
sk <- function(x) { 
  #computes the sample skewness coeff. 
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 )
}

pwr_beta = function(a){
 alpha = 0.1
 n = 20
 m = 1e4
 N = length(a)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 
 for (j in 1:N) { 
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rbeta(n, a[j], a[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr * (1-pwr) / m) 
 return(list(pwr = pwr,se = se))
}

 a = c(seq(0,1,0.1),seq(1,20,1),seq(20,100,10))
 pwr = pwr_beta(a)$pwr
 # plot the power
 se = pwr_beta(a)$se
 plot(a, pwr, type = "b", xlab = "a", ylab = "pwr", pch=16)
 abline(h = 0.1, lty = 2)
 lines(a, pwr+se, lty = 4)
 lines(a, pwr-se, lty = 4)



## ----beta---------------------------------------------------------------------

# t(v)
pwr_t = function(v){
 
 alpha = 0.1
 n = 20
 m = 1e3
 N = length(v)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 
 for (j in 1:N) { 
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rt(n,v[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr*(1-pwr) / m) 
  return(list(pwr = pwr,se = se))
}

v = seq(1,20)
pwr = pwr_t(v)$pwr
se = pwr_t(v)$se
# plot the power
plot(v, pwr, type = "b", xlab = "v", ylab = "pwr", ylim = c(0,1),pch=16)
abline(h = 0.1, lty = 2)
lines(v, pwr+se, lty = 4)
lines(v, pwr-se, lty = 4)


## -----------------------------------------------------------------------------
count5test <- function(x, y) {
        X <- x - mean(x)
        Y <- y - mean(y)
        outx <- sum(X > max(Y)) + sum(X < min(Y))
        outy <- sum(Y > max(X)) + sum(Y < min(X))
        return(as.integer(max(c(outx, outy)) > 5))
}
set.seed(12345)
alpha.hat <- 0.055
n <- c(10, 20, 50, 100, 500, 1000)
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
m <- 1e4
result <- matrix(0, length(n), 2)
for (i in 1:length(n)){
  ni <- n[i]
  tests <- replicate(m, expr={
    x <- rnorm(ni, mu1, sigma1)
    y <- rnorm(ni, mu2, sigma2)
    Fp <- var.test(x, y)$p.value
    Ftest <- as.integer(Fp <= alpha.hat)
    c(count5test(x, y), Ftest)
    })
  result[i, ] <- rowMeans(tests)
}
data.frame(n=n, C5=result[, 1], Fp=result[, 2])


## -----------------------------------------------------------------------------
library(MASS)
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  for(i in 1:c){
    central[,i]<-mydata[,i]-mean(mydata[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}

## -----------------------------------------------------------------------------
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(7912)
set.seed(7912)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=1000
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap)
n <- nrow(law)
y <- law$LSAT
z <- law$GPA
theta.hat <- cor(y,z)
theta.jack <- numeric(n)
for (i in 1:n)  theta.jack[i] <- cor(y[-i],z[-i])
bias.jack <- (n - 1) * (mean(theta.jack) - theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))

bias.jack
se.jack


## -----------------------------------------------------------------------------
library(boot) #for boot and boot.ci
data(aircondit, package = "bootstrap")
theta.boot <- function(dat, ind) {
  #function to compute the statistic
  mean(dat[ind, 1])
}
boot.obj <- boot(aircondit, statistic = theta.boot, R = 2000)
print(boot.obj)
print(boot.ci(boot.obj, type = c("norm", "basic", "perc", "bca")))


## -----------------------------------------------------------------------------
set.seed(123)
library(bootstrap)
library(boot)

lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat)
B <- 200 # number of bootstrap samples
n <- nrow(scor) # number of rows (data size)


# Jackknife
theta_j <- rep(0, n) 
for (i in 1:n) {
  x <- scor [-i,]
  lambda <- eigen(cov(x))$values
  theta_j[i] <- lambda[1] / sum(lambda)
# the i-th entry of theta_j is the i-th "leave-one-out" estimation of theta
}
bias_jack <- (n - 1) * (mean(theta_j) - theta_hat) 
# the estimated bias of theta_hat, using jackknife
se_jack <- (n - 1) * sqrt(var(theta_j) / n)
# the estimated se of theta_hat, using jackknife

# print the answers 
bias_jack
se_jack

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)


for (k in 1:n) {
  if(k < n) 
    i <- k+1  else  i <- 1
  
  y <- magnetic[c(-k,-i)]
  x <- chemical[c(-k,-i)]
 
  J1 <- lm(y ~ x)
  yhat1k <- J1$coef[1] + J1$coef[2] * chemical[k]
  yhat1i <- J1$coef[1] + J1$coef[2] * chemical[i]
  e1[k] <- (magnetic[k] + magnetic[i]- yhat1k -yhat1i)/2
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2k <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
  yhat2i <- J2$coef[1] + J2$coef[2] * chemical[i] + J2$coef[3] * chemical[i]^2
  e2[k] <- (magnetic[i] + magnetic[k] - yhat2i -yhat2k)/2
    
  J3 <- lm(log(y) ~ x)
  logyhat3k <- J3$coef[1] + J3$coef[2] * chemical[k]
  logyhat3i <- J3$coef[1] + J3$coef[2] * chemical[i]
  yhat3k <- exp(logyhat3k)
  yhat3i <- exp(logyhat3i)
  e3[k] <- (magnetic[k] + magnetic[i] - yhat3k - yhat3i)/2
    
  J4 <- lm(log(y) ~ log(x))
  logyhat4k <- J4$coef[1] + J4$coef[2] * log(chemical[k])
  logyhat4i <- J4$coef[1] + J4$coef[2] * log(chemical[i])
  yhat4k <- exp(logyhat4k)
  yhat4i <- exp(logyhat4i)
  e4[k] <- (magnetic[k] + magnetic[i] - yhat4k - yhat4i)/2
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## ---- include=FALSE-----------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
yhat4 <- exp(logyhat4)
e4[k] <- magnetic[k] - yhat4
}

## -----------------------------------------------------------------------------
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y)) 
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer((outx > 4)||(outy > 7)))
}

n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0 
sigma1 <- sigma2 <- 1

x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean 
y <- y - mean(y)
  
R <- 999
z <- c(x, y)
K <- 1:50
D <- numeric(R)

D0 <- count5test(x, y) 
for (i in 1:R) {
   #generate indices k for the first sample
  k <- sample(K, size = 20, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  D[i] <- count5test(x1, y1) 
}

mean(D)



## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(Ball)
library(energy)

m <- 1e3

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column? 
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)

alpha <- 0.1

## -----------------------------------------------------------------------------
m <- 1e3; k<-3; p<-2; mu <- 0; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)

for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
m <- 1e3; k<-3; p<-2; mu <- 0.3; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)

for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}

pow <- colMeans(p.values<alpha)
pow

