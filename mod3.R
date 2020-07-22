library("car")
data("Leinhardt")
?Leinhardt
head(Leinhardt)
str(Leinhardt)
pairs(Leinhardt)

plot(infant ~ income, data = Leinhardt)
hist(Leinhardt$income)
hist(Leinhardt$infant)
#######################
Leinhardt$loginfant = log(Leinhardt$infant)
Leinhardt$logincome = log(Leinhardt$income)
plot(loginfant ~ logincome, data = Leinhardt)
###Modeling
lmod = lm(loginfant ~ logincome, data = Leinhardt)
summary(lmod)

dat = na.omit(Leinhardt)
Leinhardt
library("rjags")
mod1_string = "model{
  for(i in 1:n){
    y[i] ~ dnorm(mu[i], prec)
    mu[i] = b[1] + b[2]*log_income[i]
  }
  for(j in 1:2){
    b[j] ~ dnorm(0.0, 1.0/1.0e6)
  }
  
  prec ~ dgamma(5.0/2.0, 5.0*10.0/2.0)
  sig2 = 1.0/prec
  sig = sqrt(sig2)
}"

set.seed(72)
data1_jags = list(y = dat$loginfant, n=nrow(dat), log_income = dat$logincome)

params1 = c("b", "sig")
inits1 = function()
{
  inits= list("b"=rnorm(2, 0.0, 100.0), "prec"=rgamma(1,1.0,1.0))
}

mod1 = jags.model(textConnection(mod1_string), data = data1_jags, inits = inits1, n.chains = 3)

update(mod1, 1000)

mod1_sim = coda.samples(model = mod1, variable.names = params1, n.iter = 5e3)

mod1_csim = do.call(rbind, mod1_sim)

###cONVERGENCE

plot(mod1_sim)

gelman.diag(mod1_sim) # Gelman and Rubin Diagnostic

autocorr.plot(mod1_sim)

effectiveSize(mod1_sim)
summary(mod1_sim)
summary(lmod)

##Residuals

lmod0 = lm(infant ~ income, data = Leinhardt) y
plot(resid(lmod0))
plot(predict(lmod0), resid(lmod0))

qqnorm(resid(lmod0))

X = cbind(rep(1.0, data1_jags$n), data1_jags$log_income)
head(X)

pm_params1 = colMeans(mod1_csim)

yhat1 = drop(X %*% pm_params1[1:2])
resid1 = data1_jags$y - yhat1
plot(resid1) # against data index
plot(yhat1, resid1) # against predicted values
qqnorm(resid1)
plot(predict(lmod), resid(lmod)) # to compare with reference linear model

rownames(dat)[order(resid1, decreasing=TRUE)[1:5]] # which countries have the largest positive residuals?
