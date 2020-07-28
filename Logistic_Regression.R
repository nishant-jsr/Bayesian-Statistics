library("boot")
data("urine")
?urine
head(urine)

dat = na.omit(urine)
dim(dat)
pairs(dat)

library("corrplot")
Cor = cor(dat)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

X = scale(dat[,-1], center=TRUE, scale=TRUE) #Variable Selection
head(X[,"gravity"])
colMeans(X)
apply(X, 2, sd)
#model

ddexp = function(x, mu, tau) {
  0.5*tau*exp(-tau*abs(x-mu)) 
}
curve(ddexp(x, mu=0.0, tau=1.0), from=-5.0, to=5.0, ylab="density", main="Double exponential\ndistribution") # double exponential distribution
curve(dnorm(x, mean=0.0, sd=1.0), from=-5.0, to=5.0, lty=2, add=TRUE) # normal distribution
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")

library("rjags")
mod1_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:6) {
        b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
    }
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"], ph=X[,"ph"], osmo=X[,"osmo"], cond=X[,"cond"], urea=X[,"urea"], calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)
summary(mod1_sim)

par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))
colnames(X) # variable names

mod2_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:3) {
        b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression
    }
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)

update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

dic2 = dic.samples(mod2, n.iter=1e3)

dic1
dic2
summary(mod2_sim)
HPDinterval(mod2_csim)

par(mfrow=c(3,1))
densplot(mod2_csim[,1:3], xlim=c(-3.0, 3.0))
colnames(X)[c(1,4,6)] # variable names

#How do we turn model parameter estimates into model predictions? The key is the form of the model.
#Remember that the likelihood is Bernoulli, which is 1 with probability p. 
#We modeled the logit of p as a linear model, which we showed in the first segment of this lesson leads to an exponential form for E(y)=p.
#Take the output from our model in the last segment. We will use the posterior means as point estimates of the parameters.

pm_coef = colMeans(mod2_csim)
pm_coef

pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3]
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)

plot(phat, jitter(dat$r))
#Suppose we choose a cutoff for these predicted probabilities. 
#If the model tells us the probability is higher than 0.5, 
#we will classify the observation as a 1 and if it is less than 0.5, we will classify it as a 0.
#That way the model classifies each data point.
#Now we can tabulate these classifications against the truth to see how well the model predicts the original data.

(tab0.5 = table(phat > 0.5, data_jags$y))
sum(diag(tab0.5)) / sum(tab0.5)
#The correct classification rate is about 76%, not too bad, but not great.

(tab0.3 = table(phat > 0.3, data_jags$y))
sum(diag(tab0.3)) / sum(tab0.3)


