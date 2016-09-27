## @knitr peak-density

x <- seq(-5, 5, len = 100)
y <- rnorm(25, 1, 0.5)

hist(y, probability = TRUE, xlim = c(-5, 5))
lines(x, dnorm(x, 1, 0.5))
lines(x, dnorm(x, 1, 2), col = 'red')

sum(dnorm(y, 1, 0.5, log = TRUE))
sum(dnorm(y, 1, 2, log = TRUE))

## @knitr occams-razor

# illustration of how hierarchical specification favors simpler models
#  if the data allow it

sigma <- 1
n <- 10
J <- 6
mu <- 2
tauSmall <- 0.1
tauLarge <- 3

### scenario 1 - little variation in theta ###
theta <- rnorm(J, mu, tauSmall)
ybar <- rnorm(J, theta, sigma/sqrt(n))

# pooled estimation with small tau (almost complete pooling)
thetaHat <- (ybar/sigma^2 + mu/tauSmall^2) / (1/sigma^2 + 1/tauSmall^2)
ll <- sum(dnorm(ybar, thetaHat, sigma/sqrt(n), log = TRUE))  # log-lik
lp <- sum(dnorm(thetaHat, mu, tauSmall, log = TRUE))             # log-prior
c(ll, lp, ll+lp)

# pooled estimation with large tau (limited pooling)
thetaHat <- (ybar/sigma^2 + mu/tauLarge^2) / (1/sigma^2 + 1/tauLarge^2)
ll <- sum(dnorm(ybar, thetaHat, sigma/sqrt(n), log = TRUE))  # log-lik
lp <- sum(dnorm(thetaHat, mu, tauLarge, log = TRUE))             # log-prior
c(ll, lp, ll+lp)

# conclusion: complexity penalty favors pooling (simpler model) here because improvement in likelihood with more complicated model is small

### scenario 2 - large variation in theta ###
theta <- rnorm(J, mu, tauLarge)
ybar <- rnorm(J, theta, sigma/sqrt(n))

# pooled estimation with small tau (almost complete pooling)
thetaHat <- (ybar/sigma^2 + mu/tauSmall^2) / (1/sigma^2 + 1/tauSmall^2)
ll <- sum(dnorm(ybar, thetaHat, sigma/sqrt(n), log = TRUE))  # log-lik
lp <- sum(dnorm(thetaHat, mu, tauSmall, log = TRUE))             # log-prior
c(ll, lp, ll+lp)

# pooled estimation with large tau (limited pooling)
thetaHat <- (ybar/sigma^2 + mu/tauLarge^2) / (1/sigma^2 + 1/tauLarge^2)
ll <- sum(dnorm(ybar, thetaHat, sigma/sqrt(n), log = TRUE))  # log-lik
lp <- sum(dnorm(thetaHat, mu, tauLarge, log = TRUE))             # log-prior
c(ll, lp, ll+lp)

# conclusion: data do not allow pooling and overwhelm the prior's preference for a simpler model


## @knitr ig-prior

# illustrating wacky behavior of the IG(eps, eps) prior
#   this prior is NOT recommended for variance components, particularly for random effects variances

par(mfrow = c(1,3),mai = c(.5,.5,.3,.1),mgp = c(2.2,0.8,0))
alpha = beta = .001
x1 = seq(0,1,len = 10000)
f1 = (1/gamma(alpha))*beta^(-alpha)*exp(-beta/x1)/(x1^(alpha+1))
plot(x1,f1,type = 'l',xlab = expression(tau^2),ylab = expression(pi(tau^2)),main = 'a.) prior on [0,1]',cex.lab = 1.5)
x2 = seq(0,.01,len = 10000)
f2 = (1/gamma(alpha))*beta^(-alpha)*exp(-beta/x2)/(x2^(alpha+1))
plot(x2,f2,type = 'l',xlab = expression(tau^2),ylab = expression(pi(tau^2)),main = 'b.) prior on [0,0.01]',cex.lab = 1.5)
x3 = seq(0,.001,len = 10000)
f3 = (1/gamma(alpha))*beta^(-alpha)*exp(-beta/x3)/(x3^(alpha+1))
plot(x3,f3,type = 'l',xlab = expression(tau^2),ylab = expression(pi(tau^2)),main = 'c.) prior on [0,0.001]',cex.lab = 1.5)

## @knitr ig-sensitivity

# some illustrations of sensitivity of posterior to prior specification with
#   IG priors for random effects variance term
#   when there is little or no variation in the random effects

# likely not demonstrated in class but might be shown in Unit 6 on computation

library(nimble)
J <- 10
n <- 10
sigma <- 1
mu <- 2
y <- matrix(rnorm(n*J, mu, sigma), nrow = n, ncol = J)
burnin <- 1000
nIts <- 100000
postBurn <- (burnin+1):nIts

code_flattau <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n[j]) {
            y[i,j] ~ dnorm(theta[j], sd = sigma)
        }
        theta[j] ~ dnorm(mu, sd = tau)
    }
    tau ~ dunif(0, 1000)
    sigma ~ dunif(0, 1000)
    mu ~ dnorm(0, sd = 1000)
})
m <- nimbleModel(code_flattau, data = list(y = y), constants = list(n = rep(n, J), J = J),
                 inits = list(mu = 1, sigma = 0.5, tau = 1))
cm <- compileNimble(m)
conf <- configureMCMC(m, monitors = c('mu', 'theta', 'sigma', 'tau'))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(nIts)
smp1 <- as.matrix(cmcmc$mvSamples)

tsplot(smp1[postBurn, 'tau'])
hist(smp1[postBurn, 'tau'], ncl = 100)
hist(smp1[postBurn, 'tau'], ncl = 100, xlim = c(0,.2))

# NIMBLE (for now) and JAGS/BUGS do not have an IG distribution
# could easily write an inverse-gamma distribution as a user-defined distribution in NIMBLE
code_igeps <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n[j]) {
            y[i,j] ~ dnorm(theta[j], sd = sigma)
        }
        theta[j] ~ dnorm(mu, itau2)
    }
    itau2 ~ dgamma(eps, eps)
    sigma ~ dunif(0, 1000)
    mu ~ dnorm(0, sd = 1000)
})

m <- nimbleModel(code_igeps, data = list(y = y), constants = list(n = rep(n, J), J = J),
                 inits = list(eps = .001, mu = 1, sigma = 0.5, itau2 = 1))
cm <- compileNimble(m)
conf <- configureMCMC(m)
#conf$removeSamplers('tau2')
#conf$addSampler('tau2', control = list(log = TRUE))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(nIts)
smp2 <- as.matrix(cmcmc$mvSamples)

tsplot(sqrt(1/smp2[postBurn, 'itau2']))
hist(sqrt(1/smp2[postBurn, 'itau2']), ncl = 100)
hist(sqrt(1/smp2[postBurn, 'itau2']), ncl = 100, xlim = c(0, .2))

pri <- sqrt(1/rgamma(10000,.001,.001))

m <- nimbleModel(code_igeps, data = list(y = y), constants = list(n = rep(n, J), J = J),
                 inits = list(eps = .00001, mu = 1, sigma = 0.5, itau2 = 1))
cm <- compileNimble(m)
conf <- configureMCMC(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(nIts)
smp3 <- as.matrix(cmcmc$mvSamples)

tsplot(sqrt(1/smp3[postBurn, 'itau2']))
hist(sqrt(1/smp3[postBurn, 'itau2']), ncl = 100)
hist(sqrt(1/smp3[postBurn, 'itau2']), ncl = 100, xlim = c(0,.2))


digimproper <- nimbleFunction(
    run = function(x = double(0), 
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- -log(x)
        if(log) return(logProb)
        else return(exp(logProb)) 
        })


rigimproper <- nimbleFunction(
    run = function(n = integer(0)) {
        returnType(double(0))
        if(n != 1) print("rigimproper only allows n = 1; using n = 1.")
        return(0)
    })

registerDistributions(list(
        digimproper = list(
               BUGSdist = "digimproper()",
               pqAvail = FALSE,
               range = c(0, Inf)
               )))

code_impropertau <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n[j]) {
            y[i,j] ~ dnorm(theta[j], sd = sigma)
        }
        theta[j] ~ dnorm(mu, var = tau2)
    }
    tau2 ~ digimproper()
    sigma ~ dunif(0, 1000)
    mu ~ dnorm(0, sd = 1000)    
})


m <- nimbleModel(code_impropertau, data = list(y = y), constants = list(n = rep(n, J), J = J),
                 inits = list(mu = 1, sigma = 0.5, tau2 = 1))
cm <- compileNimble(m)
conf <- configureMCMC(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(nIts)
smp4 <- as.matrix(cmcmc$mvSamples)

tsplot(sqrt(smp4[postBurn, 'tau2']))
hist(sqrt(smp4[postBurn, 'tau2']), ncl = 100)
hist(sqrt(smp4[postBurn, 'tau2']), ncl = 100, xlim = c(0,.2))

code_halft <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n[j]) {
            y[i,j] ~ dnorm(theta[j], sd = sigma)
        }
        theta[j] ~ dnorm(mu, sd = tau)
    }
    tau ~ T(dt(mu = 0, sigma = 1, df = 1), 0, Inf)
    sigma ~ dunif(0, 1000)
    mu ~ dnorm(0, sd = 1000)
})
m <- nimbleModel(code_halft, data = list(y = y), constants = list(n = rep(n, J), J = J),
                 inits = list(mu = 1, sigma = 0.5, tau = 1))
cm <- compileNimble(m)
conf <- configureMCMC(m, monitors = c('mu', 'theta', 'sigma', 'tau'))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(nIts)
smp5 <- as.matrix(cmcmc$mvSamples)
