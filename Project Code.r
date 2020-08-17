###ACTL30004 Assignment Group 7
library(CASdatasets)
library(maxLik)
library(rootSolve)
library(MASS)
library(countreg)
library(statmod)
library(stats)
data("ausprivauto0405")

##Part(a)
claimnb <- ausprivauto0405$ClaimNb

#Fitting a poisson distribution
lambda0 <- mean(claimnb)

loglikfun.pois <- function(param){
  lambda <- param[1]
  sum(dpois(claimnb, lambda, log=TRUE))
  }
mle.pois <- maxLik(logLik=loglikfun.pois, start=c(lambda=lambda0))
summary.pois <- summary(mle.pois)

summary.pois

###same

#Fitting a negative binomial distribution
negbin <- function(y,r,beta){
  factorial(y+r-1)/(factorial(y)*factorial(r-1))*(1/(1+beta))^r*(beta/(1+beta))^y
  }
###same

beta0 <- var(claimnb)/mean(claimnb) - 1
r0 <- mean(claimnb)/beta0
### initial est of mme
### estimates are still same

loglikfun.negbin <- function(param){
  r <- param[1]
  beta <- param[2]
  sum(log(negbin(claimnb, r, beta)))
}
mle.negbin <- maxLik(logLik=loglikfun.negbin, start=c(r=r0, beta=beta0))
summary.negbin <- summary(mle.negbin)

summary.negbin

### mle estimates vary slightly based on method of initial starting values??? 

#Fitting a zero inflated poisson distribution
mme.zipois <- function(x){
  E <-(1-x[2])*x[1]-mean(claimnb)
  Var <- x[1]*(1-x[2])*(1+x[2]*x[1])-var(claimnb)
  c(E,Var)
}

initial.zipois <- multiroot(mme.zipois, start=c(0.1,0.5)) 
theta0 <- initial.zipois$root[1]
psi0 <- initial.zipois$root[2]

###same

loglikfun.zipois <- function(param){
  theta <- param[1]
  psi <- param[2]
  sum(dzipois(claimnb, theta, psi, log=TRUE))
  }
mle.zipois <- maxLik(logLik=loglikfun.zipois, start=c(theta=theta0,psi=psi0))
summary.zipois <- summary(mle.zipois)

summary.zipois

loglik.pois <- summary.pois$loglik
loglik.negbin <- summary.negbin$loglik
loglik.zipois <- summary.zipois$loglik

#Likelihood Ratio Test (LRT) @ 5% significance level
#Poisson vs NB
tstat1 <- 2*(loglik.negbin-loglik.pois)
ifelse(tstat1 > qchisq(0.05, 1, lower.tail=FALSE), "Negative binomial model is preferable", "Poisson model is preferable")
  
#Poisson vs ZIPois
tstat2 <- 2*(loglik.zipois-loglik.pois)
ifelse(tstat2 > qchisq(0.05, 1, lower.tail=FALSE), "Zero inflated poisson model is preferable", "Poisson model is preferable")

#can't do LRT for NB vs ZIPois since not nested

#Akaike's Information Criterion (AIC)
AIC.pois <- -2*loglik.pois + 2*1
AIC.negbin <- -2*loglik.negbin + 2*2
AIC.zipois <- -2*loglik.zipois + 2*2
if(AIC.pois < AIC.negbin && AIC.pois < AIC.zipois){
  print("Poisson model is preferable")
  } else if(AIC.negbin < AIC.zipois){
  print("Negative binomial model is preferable")
  } else{
  print("Zero inflated poisson model is preferable")
  }

#Bayesian Information Criterion (BIC)
BIC.pois <- -2*loglik.pois + 1*log(67856)
BIC.negbin <- -2*loglik.negbin + 2*log(67856)
BIC.zipois <- -2*loglik.zipois + 2*log(67856)
if(BIC.pois < BIC.negbin && BIC.pois < BIC.zipois){
  print("Poisson model is preferable")
  } else if(AIC.negbin < AIC.zipois){
  print("Negative binomial model is preferable")
  } else{
  print("Zero inflated poisson model is preferable")
  }

###same

##Part(b)

##maybe use attach if can be bothered
vehvalue <- ausprivauto0405$VehValue
vehage <- as.integer(ausprivauto0405$VehAge %in% c("old cars", "oldest cars"))
drivage <- as.integer(ausprivauto0405$DrivAge %in% c("old people", "older work. people", "oldest people"))
exposure <- ausprivauto0405$Exposure

#Fitting a poisson GLM
#for starting values
poislm <- lm(claimnb ~ vehvalue + vehage + drivage, offset=log(exposure))

poislm

poisglm <- glm(claimnb ~ vehvalue + vehage + drivage, offset=log(exposure), family=poisson(link=log), start=coef(poislm))
summary.poisglm <- summary(poisglm)

summary.poisglm

###same

##Part(e)
intercept <- c(rep(1, 67856))
negbin.repar <- function(y,r,mu){
  factorial(y + r - 1)/(factorial(y)*factorial(r - 1))*(r/(r + mu))^r*(mu/(r + mu))^y
}

#Initial values
new.r <- coef(summary.negbin)[1, 1]
new.a <- coef(summary.poisglm)[1, 1]
new.b <- coef(summary.poisglm)[2, 1]
new.c <- coef(summary.poisglm)[3, 1]
new.d <- coef(summary.poisglm)[4, 1]
new.mu <- exp(new.a + new.b*vehvalue + new.c*vehage + new.d*drivage + log(exposure))

#Fisher-scoring algorithm
U1 <- 1
U2 <- 1 
U3 <- 1
U4 <- 1
i <- 1
tol.level <- 1e-10

while((abs(U1) > tol.level) | (abs(U2) > tol.level) | (abs(U3) > tol.level) | (abs(U4) > tol.level)){
  r0 <- new.r
  a0 <- new.a
  b0 <- new.b
  c0 <- new.c
  d0 <- new.d
  mu0 <- new.mu
  
  I11 <- sum(r0*mu0*intercept*intercept/(r0 + mu0))
  I12 <- sum(r0*mu0*intercept*vehvalue/(r0 + mu0))
  I13 <- sum(r0*mu0*intercept*vehage/(r0 + mu0))
  I14 <- sum(r0*mu0*intercept*drivage/(r0 + mu0))
  I22 <- sum(r0*mu0*vehvalue*vehvalue/(r0 + mu0))
  I23 <- sum(r0*mu0*vehvalue*vehage/(r0 + mu0))
  I24 <- sum(r0*mu0*vehvalue*drivage/(r0 + mu0))
  I33 <- sum(r0*mu0*vehage*vehage/(r0 + mu0))
  I34 <- sum(r0*mu0*vehage*drivage/(r0 + mu0))
  I44 <- sum(r0*mu0*drivage*drivage/(r0 + mu0))
  fim <- matrix(c(I11, I12, I13, I14, I12, I22, I23, I24, I13, I23, I33, I34, I14, I24, I34, I44), 4, 4)
  inv.fim=solve(fim)
  
  U1 <- sum((r0*(claimnb - mu0))*intercept/(r0 + mu0))
  U2 <- sum((r0*(claimnb - mu0))*vehvalue/(r0 + mu0))
  U3 <- sum((r0*(claimnb - mu0))*vehage/(r0 + mu0))                
  U4 <- sum((r0*(claimnb - mu0))*drivage/(r0 + mu0))   
  
  a1 <- a0 + inv.fim[1, 1]*U1 + inv.fim[1, 2]*U2 + inv.fim[1, 3]*U3 + inv.fim[1, 4]*U4
  b1 <- b0 + inv.fim[2, 1]*U1 + inv.fim[2, 2]*U2 + inv.fim[2, 3]*U3 + inv.fim[2, 4]*U4
  c1 <- c0 + inv.fim[3, 1]*U1 + inv.fim[3, 2]*U2 + inv.fim[3, 3]*U3 + inv.fim[3, 4]*U4
  d1 <- d0 + inv.fim[4, 1]*U1 + inv.fim[4, 2]*U2 + inv.fim[4, 3]*U3 + inv.fim[4, 4]*U4
  mu1 <- exp(a1 + b1*vehvalue + c1*vehage + d1*drivage + log(exposure))
  
  loglikfun.r <- function(param){
    r1 <- param[1]
    sum(log(negbin.repar(claimnb, r1, mu1)))
    }
  mle.r <- maxLik(logLik=loglikfun.r, start=c(r1=r0))
  summary.negbinglm1 <- summary(mle.r)
  
  new.mu <- mu1
  scovec <- c(U1, U2, U3, U4)
  
  new.r <- coef(summary.negbinglm1)[1, 1]
  new.a <- a1
  new.b <-b1
  new.c <-c1
  new.d <- d1
  sd.r <- coef(summary.negbinglm1)[1, 2]
  sd.a <- sqrt(inv.fim[1, 1])
  sd.b <- sqrt(inv.fim[2, 2])
  sd.c <- sqrt(inv.fim[3, 3])
  sd.d <- sqrt(inv.fim[4, 4])
  
  print(paste0("Iteration number:", i))
  print("Score vector:")
  print(scovec)
  print("Maximum likelihood estimates:")
  print(paste0("r", i, "=", new.r))
  print(paste0("a", i, "=", new.a))
  print(paste0("b", i, "=", new.b))
  print(paste0("c", i, "=", new.c))
  print(paste0("d", i, "=", new.d))
  print("Standard errors:")
  print(paste0("r", i, "=", sd.r))
  print(paste0("a", i, "=", sd.a))
  print(paste0("b", i, "=", sd.b))
  print(paste0("c", i, "=", sd.c))
  print(paste0("d", i, "=", sd.d))
  print("Maximum log-likelihood:")
  print(summary.negbinglm1$loglik)
  i <- i + 1}


### estimates are very similar to original (only differ after around 5 or so dec places)
### this has 12 iterations instead of original 10

##Part(f)
loglik.m1 <- summary.negbinglm1$loglik
m2 <- glm.nb(claimnb ~ vehvalue + vehage + offset(log(exposure)))
loglik.m2 <- m2$twologlik/2
tstat3 <- 2*(loglik.m1-loglik.m2)
ifelse(tstat3 > qchisq(0.05, 1, lower.tail=FALSE), "Model is better with indicator variable for the age of the driver", "Model is not better with indicator variable for the age of the driver")

### looks good, other than changed summary.negbinglm to summary.negbinglm1


##Part(g)

#Likelihood Ratio Test (LRT) @ 5% significance level
loglik.poisglm <- -(AIC(poisglm) - 2*1)/2
loglik.negbinglm <- loglik.m1
tstat4 = 2*(loglik.negbinglm - loglik.poisglm)
ifelse(tstat4 > qchisq(0.05, 1, lower.tail=FALSE), "Negative binomial GLM is preferable", "Poisson GLM is preferable")

###looks good

#Akaike's Information Criterion (AIC)
AIC.poisglm <- AIC(poisglm)
AIC.negbinglm <- -2*loglik.negbinglm + 2*5
#changed this to loglik.negbinglm instead of indexing again
ifelse(AIC.poisglm < AIC.negbinglm, "Poisson GLM is preferable", "Negative binomial GLM is preferable")

##JUST FOR CHECKS
#negbinglm <- glm.nb(claimnb ~ vehvalue + vehage + drivage + offset(log(exposure)))
#summary.negbinglm2 <- summary(negbinglm)

### all g 

##Part(h)
mu.pois <- exp(coef(poisglm)[1] + coef(poisglm)[2]*vehvalue + coef(poisglm)[3]*vehage + coef(poisglm)[4]*drivage + log(exposure))
pearsonres.pois <- (claimnb - mu.pois)/sqrt(mu.pois)
devianceres.pois <- sign(claimnb - mu.pois)*sqrt(2*(mu.pois - claimnb*(1 - log(claimnb/mu.pois))))

#changed this to new.mu, already defined in the last iteration of e
mu.negbin <- new.mu
pearsonres.negbin <- (claimnb - mu.negbin)/sqrt(mu.negbin*(new.r + mu.negbin)/new.r)
devianceres.negbin <- sign(claimnb - mu.negbin)*sqrt(2*(new.r*log((new.r + mu.negbin)/(new.r + claimnb)) + claimnb*log((claimnb*(new.r + mu.negbin))/(mu.negbin*(new.r + claimnb)))))


par(mfrow=c(1,1))
#QQ plot for pearson's residuals of poisson GLM
qqnorm(pearsonres.pois, main="Poisson Pearson's Residuals Q-Q Plot", xlab="Standard Normal Quantiles", ylab="Residuals")
qqline(pearsonres.pois, col=2)

#QQ plot for deviance residuals of poisson GLM
qqnorm(devianceres.pois, main="Poisson Deviance Residuals Q-Q Plot", xlab="Standard Normal Quantiles", ylab="Residuals")
qqline(devianceres.pois, col=2)

#QQ plot for pearson's residuals of negative binomial GLM
qqnorm(pearsonres.negbin, main="Negative Binomial Pearson's Residuals Q-Q Plot", xlab="Standard Normal Quantiles", ylab="Residuals")
qqline(pearsonres.negbin, col=2)

#QQ plot for deviance residuals of negative binomial GLM
qqnorm(devianceres.negbin, main="Negative Binomial Deviance Residuals Q-Q Plot", xlab="Standard Normal Quantiles", ylab="Residuals")
qqline(devianceres.negbin, col=2)

### seems good

##Part(j)
id <- 1:nrow(ausprivauto0405)
ausprivauto0405$Id <- id
ausprivauto0405$Intercept <- intercept
ausprivauto0405$vehAge <- vehage
ausprivauto0405$drivAge <- drivage

#Creating datasets
training <- subset.data.frame(ausprivauto0405, (Id%%2==1), select=c(ClaimNb, Exposure, Intercept, VehValue, vehAge, drivAge))
validation <- subset.data.frame(ausprivauto0405, !(Id%%2), select=c(ClaimNb, Exposure, Intercept, VehValue, vehAge, drivAge))

#Fitting the training dataset
model.a <- glm.nb(ClaimNb ~ VehValue + vehAge + drivAge + offset(log(Exposure)), data=training)
summary.a <- summary(model.a)
summary.a

#Plotting randomized quantile residuals for validation dataset
model.b <- glm.nb(ClaimNb ~ VehValue + vehAge + drivAge + offset(log(Exposure)), data=validation)
quantileres.b <- qres.nbinom(model.b)

par(mfrow=c(1,2))
qqnorm(quantileres.b, main="Dataset B Randomised Quantile Residuals Q-Q Plot", xlab="Standard Normal Quantiles", ylab="Residuals")
qqline(quantileres.b, col=2)

#Plotting randomized quantile residuals for complete dataset
model.full <- glm.nb(claimnb ~ vehvalue + vehage + drivage + offset(log(exposure)))
quantileres.full <- qres.nbinom(model.full)

qqnorm(quantileres.full, main="Complete Dataset Randomised Quantile Residuals Q-Q Plot", xlab="Standard Normal Quantiles", ylab="Residuals")
qqline(quantileres.full, col=2)


##Part(k)
claim.amounts <- ausprivauto0405$ClaimAmount
claims.subset <- claim.amounts[claim.amounts>0]
claims.subset.small <- claims.subset[claims.subset<50000]
par(mfrow = c(1,1))

#assuming he wants the histogram of densities bc he said to superimpose densities afterwards
#however, empirical distribution implies frequencies????
hist(claims.subset.small, freq = FALSE, main = "Distribution of claims less than 50,000", xlab = "Claim Amount", ylab = "Density", col = "plum", xlim = c(0, 50000), ylim = c(0,0.0008), breaks = 80)

#lognormal
#we have closed form expression for mle for lognormal. will compute both ways to confirm answer:

#1: newton-raphson
loglikfun.lnorm <- function(param) {
  mu = param[1]
  sig = param[2]
  sum(dlnorm(claims.subset.small, meanlog = mu, sdlog = sig, log = TRUE))
}

#obtain mme for initial values
mme.lognormal = vector()
mme.lognormal[1] = 2*log(mean(claims.subset.small)) - 0.5*log(var(claims.subset.small)+(mean(claims.subset.small))^2)
mme.lognormal[2] = sqrt(log(1 + var(claims.subset.small)/((mean(claims.subset.small))^2)))

mle.lognormal=maxLik(logLik=loglikfun.lnorm,start=mme.lognormal)

#2: closed form mle
#lnorm.mu.mle = sum(log(claims.subset.small))/length(claims.subset.small)
#lnorm.mu.mle
#sqrt((sum((log(claims.subset.small) - lnorm.mu.mle)^2))/length(claims.subset.small))

#both give same values, and are close to mme values

#superimpose density on histogram
lines(0:50000, dlnorm(0:50000, meanlog = coef(mle.lognormal)[1], sdlog = coef(mle.lognormal)[2]), lwd = 2, col = "blue")

#inverse gaussian

#1: newton raphson

loglikfun.invgauss <- function(param) {
  mu.0 = param[1]
  lambda.0 = param[2]
  sum(dinvgauss(claims.subset.small, mean = mu.0, shape = lambda.0, log = TRUE))
}

mme.invgauss = vector()
mme.invgauss[1] = mean(claims.subset.small)
mme.invgauss[2] = (mean(claims.subset.small))^3/var(claims.subset.small)

mle.invgauss=maxLik(logLik=loglikfun.invgauss,start=mme.invgauss)

#closed form mle
#mean(claims.subset.small)
#length(claims.subset.small)/sum(1/claims.subset.small - 1/mean(claims.subset.small))

#both have same values for mu, lambda slightly varies with maxLik depending on initial estimates
#not sure if something is wrong with code, but would rather to closed form mle
#as it is more robust and expression is simple

#stick with closed form for now:
mle.invgauss.closed = vector()
mle.invgauss.closed[1] = mean(claims.subset.small)
mle.invgauss.closed[2] = length(claims.subset.small)/sum(1/claims.subset.small - 1/mean(claims.subset.small))

lines(0:50000, dinvgauss(0:50000, mean = mle.invgauss.closed[1], shape = mle.invgauss.closed[2]), lwd = 2, col = "red", lty = 3)
legend("topright", c("Lognormal", "Inverse Gaussian"), fill=c("blue", "red"))

#don't know what we're supposed to say for LRT

##Part(l)

#refit data including last data point
#double check that this is what he wants though

loglikfun.lnorm <- function(param) {
  mu = param[1]
  sig = param[2]
  sum(dlnorm(claims.subset, meanlog = mu, sdlog = sig, log = TRUE))
}

mme.lognormal = vector()
mme.lognormal[1] = 2*log(mean(claims.subset)) - 0.5*log(var(claims.subset)+(mean(claims.subset))^2)
mme.lognormal[2] = sqrt(log(1 + var(claims.subset)/((mean(claims.subset))^2)))

mle.lognormal=maxLik(logLik=loglikfun.lnorm,start=mme.lognormal)

mle.invgauss.closed = vector()
mle.invgauss.closed[1] = mean(claims.subset)
mle.invgauss.closed[2] = length(claims.subset)/sum(1/claims.subset - 1/mean(claims.subset))

#calculate VAR using each distribution
security.levels = c(0.9, 0.95, 0.99)

#Lognormal
lnorm.var <- qlnorm(security.levels, meanlog = coef(mle.lognormal)[1], sdlog = coef(mle.lognormal)[2])

#Inv Gaussian
invgauss.var <- qinvgauss(security.levels, mean = mle.invgauss.closed[1], shape = mle.invgauss.closed[2])

#Empirical Cdf
obs.var <- quantile(claims.subset, security.levels)

VAR <- matrix(c(lnorm.var, invgauss.var, obs.var), ncol = 3, byrow = TRUE)
colnames(VAR) <- c('0.90', '0.95', '0.99')
rownames(VAR) <- c('Lognormal', 'Inverse Gaussian', 'Observed')
VAR

#lognormal has lower VARs at each security level 
#lognormal assumes less capital required at each security level
#implies inv gaussian model is slightly more conservative
#neither of them are conservative enough for empirical

##Part(m)

#using the updated distribution to include 4624 data points

#Lognormal

#initialise
ks.lnorm = 0
sorted.claims = sort(claims.subset)
n = length(sorted.claims)

#CDF at y[i]th data point from claims distribution
F.claims <- function(i) {
  plnorm(sorted.claims[i], meanlog = coef(mle.lognormal)[1], sdlog = coef(mle.lognormal)[2])
}

#get observed statistic
for (i in 2:n) {
  ks = max(F.claims(i)-(i-1)/n, i/n-F.claims(i))
  if (ks > ks.lnorm) {
    ks.lnorm = ks
  }
}

ITER = 10000
ks.mat.lnorm = rep(0, ITER)

#CDF at y[i]th data point from simulated distribution
F.data <- function(j) {
  plnorm(sorted.data.lnorm[j], meanlog = coef(mle.lognormal)[1], sdlog = coef(mle.lognormal)[2])
}

#simulate statistics
#this takes up to 10 minutes to load on a mac, p-val will be zero tho bc maximum deviance
#is around 0.03, nowhere close to 0.1
for(i in 1:ITER) {
  #simulate and sort data
  sorted.data.lnorm = sort(rlnorm(n, meanlog = coef(mle.lognormal)[1], sdlog = coef(mle.lognormal)[2]))
  
  #calculate ks statistic for each simulation, store in matrix
  for (j in 2:n) {
    ks.lnorm = max(F.data(j)-(j-1)/n, j/n-F.data(j))
    if (ks.lnorm > ks.mat.lnorm[i]) {
      ks.mat.lnorm[i] = ks.lnorm
    }
  }
}

#simulated p val
pval.ks.lnorm <- sum(ks.mat.lnorm>=ks.lnorm)/ITER


#Inverse Gaussian

#same steps, replace lnorm with invgauss, and F with G
ks.invgauss = 0

G.claims <- function(i) {
  pinvgauss(sorted.claims[i], mean = mle.invgauss.closed[1], shape = mle.invgauss.closed[2])
}

#get observed statistic
for (i in 2:n) {
  ks = max(G.claims(i)-(i-1)/n, i/n-G.claims(i))
  if (ks > ks.invgauss) {
    ks.invgauss = ks
  }
}

ks.mat.invgauss = rep(0, ITER)

G.data <- function(j) {
  pinvgauss(sorted.data.invgauss[j], mean = mle.invgauss.closed[1], shape = mle.invgauss.closed[2])
}

for(i in 1:ITER) {
  sorted.data.invgauss = sort(rinvgauss(n, mean = mle.invgauss.closed[1], shape = mle.invgauss.closed[2]))
  
  for (j in 2:n) {
    ks.invgauss = max(G.data(j)-(j-1)/n, j/n-G.data(j))
    if (ks.invgauss > ks.mat.invgauss[i]) {
      ks.mat.invgauss[i] = ks.invgauss
    }
  }
}

pval.ks.invgauss <- sum(ks.mat.invgauss>=ks.invgauss)/ITER

#maybe we should explain that the p-val is not exactly zero, nonetheless we still need to reject our
#null hypothesis

#CHECK:
#kstest.lnorm <- ks.test(claims.subset,"plnorm", meanlog = coef(mle.lognormal)[1], sdlog = coef(mle.lognormal)[2])
#kstest.invgauss <- ks.test(claims.subset,"pinvgauss", mean = mle.invgauss.closed[1], shape = mle.invgauss.closed[2])

#observed test statistics are correct by checking kstest

#conclusion: claim amounts are neither lognormally distributed nor inverse gaussian


##Part (n)

#remember to write adv and disadv

#Lognormal
#similar process as before
#these take over 20 minutes to run
#similar to before, max value of these is around 3, observed values are around 8

ad.lnorm = n
for (i in 1:n) {
  ad.lnorm = ad.lnorm + (2*i-1)/n*(log(F.claims(i))+log(1-F.claims(n+1-i)))
}
ad.lnorm = -ad.lnorm
ad.lnorm = sqrt(ad.lnorm)

ad.mat.lnorm = rep(0, ITER)

for(i in 1:ITER) {
  sorted.data.lnorm = sort(rlnorm(n, meanlog = coef(mle.lognormal)[1], sdlog = coef(mle.lognormal)[2]))
  
  ad.mat.lnorm[i] = n
  for (j in 1:n) {
    ad.mat.lnorm[i] = ad.mat.lnorm[i]+ (2*j-1)/n*(log(F.data(j))+log(1-F.data(n+1-j)))
  }
  ad.mat.lnorm[i] = -ad.mat.lnorm[i]
  ad.mat.lnorm[i] = sqrt(ad.mat.lnorm[i])
}

pval.ad.lnorm <- sum(ad.mat.lnorm>=ad.lnorm)/ITER

#from examples i've seen online, the simulated and actual AD statistics do seem reasonable
#idk how exactly to check though

#Inverse Gaussian 

ad.invgauss = n
for (i in 1:n) {
  ad.invgauss = ad.invgauss + (2*i-1)/n*(log(G.claims(i))+log(1-G.claims(n+1-i)))
}
ad.invgauss = -ad.invgauss
ad.invgauss = sqrt(ad.invgauss)

ad.mat.invgauss = rep(0, ITER)

for(i in 1:ITER) {
  sorted.data.invgauss = sort(rinvgauss(n, mean = mle.invgauss.closed[1], shape = mle.invgauss.closed[2]))
  
  ad.mat.invgauss[i] = n
  for (j in 1:n) {
    ad.mat.invgauss[i] = ad.mat.invgauss[i]+ (2*j-1)/n*(log(G.data(j))+log(1-G.data(n+1-j)))
  }
  ad.mat.invgauss[i] = -ad.mat.invgauss[i]
  ad.mat.invgauss[i] = sqrt(ad.mat.invgauss[i])
}

pval.ad.invgauss <- sum(ad.mat.invgauss>=ad.invgauss)/ITER


