#Let’s explore an example with two factors. We’ll use the Warpbreaks data set in R. 
#Check the documentation for a description of the data by typing ?warpbreaks.

data("warpbreaks")
?warpbreaks
head(warpbreaks)

table(warpbreaks$wool, warpbreaks$tension)


boxplot(breaks ~ wool + tension, data=warpbreaks)


boxplot(log(breaks) ~ wool + tension, data=warpbreaks)


#The different groups have more similar variance if we use the logarithm of breaks. 
#From this visualization, it looks like both factors may play a role in the number of breaks. 
#It appears that there is a general decrease in breaks as we move from low to medium to high tension. 
#Let’s start with a one-way model using tension only.

#One-way model
library("rjags")

mod1_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[tensGrp[i]], prec)
    }
    
    for (j in 1:3) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(5/2.0, 5*2.0/2.0)
    sig = sqrt(1.0 / prec)
} "

set.seed(83)
str(warpbreaks)

data1_jags = list(y=log(warpbreaks$breaks), tensGrp=as.numeric(warpbreaks$tension))

params1 = c("mu", "sig")

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5e3)

## convergence diagnostics
plot(mod1_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)
Here are the results.

summary(mod1_sim)

#The 95% posterior interval for the mean of group 2 (medium tension) overlaps with both the low and high groups,
# but the intervals for low and high group only slightly overlap. 
#That is a pretty strong indication that the means for low and high tension are different. 
#Let’s collect the DIC for this model and move on to the two-way model.

dic1 = dic.samples(mod1, n.iter=1e3)

#Two-way additive model
#With two factors, one with two levels and the other with three, we have six treatment groups, which is the same situation we discussed when introducing multiple factor ANOVA. 
#We will first fit the additive model which treats the two factors separately with no interaction. 
#To get the X matrix (or design matrix) for this model, we can create it in R.

X = model.matrix( ~ wool + tension, data=warpbreaks)
head(X)

tail(X)

#By default, R has chosen the mean for wool A and low tension to be the intercept. 
#Then, there is an effect for wool B, and effects for medium tension and high tension, 
#each associated with dummy indicator variables.

mod2_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = int + alpha*isWoolB[i] + beta[1]*isTensionM[i] + beta[2]*isTensionH[i]
    }
    
    int ~ dnorm(0.0, 1.0/1.0e6)
    alpha ~ dnorm(0.0, 1.0/1.0e6)
    for (j in 1:2) {
        beta[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(3/2.0, 3*1.0/2.0)
    sig = sqrt(1.0 / prec)
} "

data2_jags = list(y=log(warpbreaks$breaks), isWoolB=X[,"woolB"], isTensionM=X[,"tensionM"], isTensionH=X[,"tensionH"])

params2 = c("int", "alpha", "beta", "sig")

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

## convergence diagnostics
plot(mod2_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)
#Let’s summarize the results, collect the DIC for this model, and compare it to the first one-way model.

summary(mod2_sim)

(dic2 = dic.samples(mod2, n.iter=1e3))

dic1

#This suggests there is much to be gained adding the wool factor to the model. 
#Before we settle on this model however, we should consider whether there is an interaction. 
#Let’s look again at the box plot with all six treatment groups.

boxplot(log(breaks) ~ wool + tension, data=warpbreaks)


#Our two-way model has a single effect for wool B and the estimate is negative. 
#If this is true, then we would expect wool B to be associated with fewer breaks than its wool A counterpart on average. 
#This is true for low and high tension, but it appears that breaks are higher for wool B when there is medium tension. 
#That is, the effect for wool B is not consistent across tension levels,
# so it may appropriate to add an interaction term. In R, this would look like:
  
  lmod2 = lm(log(breaks) ~ .^2, data=warpbreaks)
summary(lmod2)

#Adding the interaction, we get an effect for being in wool B and medium tension, 
#as well as for being in wool B and high tension.
# There are now six parameters for the mean, one for each treatment group,
# so this model is equivalent to the full cell means model. Let’s use that.

#Two-way cell means model
#In this new model, μ will be a matrix with six entries, each corresponding to a treatment group.

mod3_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[woolGrp[i], tensGrp[i]], prec)
    }
    
    for (j in 1:max(woolGrp)) {
        for (k in 1:max(tensGrp)) {
            mu[j,k] ~ dnorm(0.0, 1.0/1.0e6)
        }
    }
    
    prec ~ dgamma(3/2.0, 3*1.0/2.0)
    sig = sqrt(1.0 / prec)
} "

str(warpbreaks)

data3_jags = list(y=log(warpbreaks$breaks), woolGrp=as.numeric(warpbreaks$wool), tensGrp=as.numeric(warpbreaks$tension))

params3 = c("mu", "sig")

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, n.chains=3)
update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

plot(mod3_sim, ask=TRUE)

## convergence diagnostics
gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
effectiveSize(mod3_sim)
raftery.diag(mod3_sim)
#Let’s compute the DIC and compare with our previous models.

(dic3 = dic.samples(mod3, n.iter=1e3))

dic2

dic1

#This suggests that the full model with interaction between wool and tension (which is equivalent to the cell means model)
# is the best for explaining/predicting warp breaks.

Results
summary(mod3_sim)

HPDinterval(mod3_csim)

par(mfrow=c(3,2)) # arrange frame for plots
densplot(mod3_csim[,1:6], xlim=c(2.0, 4.5))


##It might be tempting to look at comparisons between each combination of treatments, but we warn that this could yield spurious results. 
##When we discussed the statistical modeling cycle, we said it is best not to search your results for interesting hypotheses, 
#because if there are many hypotheses, some will appear to show “effects” or “associations” simply due to chance. 
##Results are most reliable when we determine a relatively small number of hypotheses we are interested in beforehand, collect the data, and statistically evaluate the evidence for them.

#One question we might be interested in with these data is finding the treatment combination that produces the fewest breaks. 
#To calculate this, we can go through our posterior samples and for each sample, 
#find out which group has the smallest mean. 
#These counts help us determine the posterior probability that each of the treatment groups has the smallest mean.

prop.table( table( apply(mod3_csim[,1:6], 1, which.min) ) )
## 
##          2          3          4          5          6 
## 0.01553333 0.11880000 0.01060000 0.11233333 0.74273333
#The evidence supports wool B with high tension as the treatment that produces the fewest breaks.