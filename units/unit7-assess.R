## @knitr schools-hyperprior

library(nimble)
library(methods) # some issue with using Nimble in contexts in which methods is not loaded

y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)
J <- length(y)

## "non-informative" IG priors

## assess sensitivity to different IG priors

codeIG <- nimbleCode({
    for(j in 1:J) {
        y[j] ~ dnorm(theta[j], sd = sigma[j])
        theta[j] ~ dnorm(mu, itau2)
    }
    mu ~ dnorm(0, .000001)
    itau2 ~ dgamma(eps, eps)
    tau <- 1/sqrt(itau2)
})

nIts <- 5000
burnin <- 1000

set.seed(0)
eps <- c(1, .001)

modelIG <- cmodelIG <- mcmcIG <- cmcmcIG <- smpIG <- list(NULL, NULL)
for(prior in 1:2) {
    modelIG[[prior]] <- nimbleModel(codeIG, data = list(y = y),
                                constants = list(eps = eps[prior],
                                                 sigma = sigma, J = J),
                                inits = list(theta = y, mu = mean(y),
                                             itau2 = 1/var(y)))

    cmodelIG[[prior]] <- compileNimble(modelIG[[prior]])
    mcmcIG[[prior]] <- buildMCMC(modelIG[[prior]], monitors = c('theta', 'mu', 'tau'))
    cmcmcIG[[prior]] <- compileNimble(mcmcIG[[prior]], project = modelIG[[prior]])
    cmcmcIG[[prior]]$run(nIts)
    smpIG[[prior]] <- as.matrix(cmcmcIG[[prior]]$mvSamples)[(burnin+1):nIts, ]
}

tau <- seq(0, 30, length = 100000)
# IG density for tau^2 translated to density for tau
dtau <- function(tau, eps) {
    alpha <- beta <- eps
    return(exp(log(2) + log(tau) + alpha*log(beta) - lgamma(alpha) -(alpha+1)*log(tau^2) - beta / tau^2))
}


par(mfrow = c(2, 2))
hist(smpIG[[1]][ , 'tau'], xlim = c(0, 30), breaks = 30, main = "IG(1,1)",
     xlab = expression(tau), probability = TRUE)
lines(tau, sapply(tau, dtau, 1))
hist(smpIG[[2]][ , 'tau'], xlim = c(0, 30), breaks = 60, main = "IG(.001, .001)",
     xlab = expression(tau), probability = TRUE)
# why is BDA Fig. 5.9a showing different result? also, my density doesn't appear to integrate to 1... 
lines(tau, sapply(tau, dtau, .001))
ts.plot(smpIG[[1]][ , 'tau'], main = "IG(1,1)",
     ylab = expression(tau))
ts.plot(smpIG[[2]][ , 'tau'], main = "IG(.001, .001)",
     ylab = expression(tau))

## state-of-the-art non-informative priors

## assess sensitivity to uniform and half-t priors for tau
codeGelman <- nimbleCode({
    for(j in 1:J) {
        y[j] ~ dnorm(theta[j], sd = sigma[j])
        theta[j] ~ dnorm(mu, sd = tau)
    }
    mu ~ dnorm(0, .000001)
    if(prior == 1)
        tau ~ dunif(0, 50)
    else tau ~ T(dt(0, sigma = 25, df = 1), 0, Inf)
})

nIts <- 50000
set.seed(0)

modelGelman <- cmodelGelman <- confGelman <- mcmcGelman <- cmcmcGelman <- smpGelman <- list(NULL, NULL)
for(prior in 1:2) {
    modelGelman[[prior]] <- nimbleModel(codeGelman, data = list(y = y),
                                constants = list(sigma = sigma, J = J),
                                inits = list(theta = y, mu = mean(y),
                                             tau = sd(y)))

    cmodelGelman[[prior]] <- compileNimble(modelGelman[[prior]])
    confGelman[[prior]] <- configureMCMC(modelGelman[[prior]], monitors = c('theta', 'mu', 'tau'))
    confGelman[[prior]]$removeSamplers('tau')
    confGelman[[prior]]$addSampler('tau', type = 'RW', control = list(log = TRUE))
    mcmcGelman[[prior]] <- buildMCMC(confGelman[[prior]])
    cmcmcGelman[[prior]] <- compileNimble(mcmcGelman[[prior]], project = modelGelman[[prior]])
    cmcmcGelman[[prior]]$run(nIts)
    smpGelman[[prior]] <- as.matrix(cmcmcGelman[[prior]]$mvSamples)[(burnin+1):nIts, ]
}

par(mfrow = c(3, 2))
hist(smpGelman[[1]][ , 'tau'], xlim = c(0, 30), breaks = 30, main = "uniform",
     xlab = expression(tau))
hist(smpGelman[[2]][ , 'tau'], xlim = c(0, 30), breaks = 30, main = "half-t",
     xlab = expression(tau))
ts.plot(smpGelman[[1]][ , 'tau'], main = "uniform",
     ylab = expression(tau), ylim = c(0, 50))
ts.plot(smpGelman[[2]][ , 'tau'], main = "half-t",
     ylab = expression(tau), ylim = c(0, 50))
# sensitivity of a parameter we might specifically care about?
hist(smpGelman[[1]][ , 'mu'], xlim = c(0, 30), breaks = 60,
     main = "mu posterior with uniform on tau",
     xlab = expression(mu))
hist(smpGelman[[2]][ , 'mu'], xlim = c(0, 30), breaks = 30, main = "mu posterior with half-t on tau",
     xlab = expression(mu))

quantile(smpGelman[[1]][ , 'mu'], c(.025, .5, .975))
quantile(smpGelman[[2]][ , 'mu'], c(.025, .5, .975))



## @knitr schools-hyperprior-analytic

# check the result for what mode is using our gridded marginal for tau
dtau <- function(tau, y, sigma, unif = TRUE) {
    vInvHat <- sum(1 / (sigma^2 + tau^2))
    muHat <- sum(y/(sigma^2 + tau^2)) / vInvHat
    terms <- -0.5*log(vInvHat) -0.5*sum(log(sigma^2+tau^2)) -0.5*sum((y - muHat)^2 / (sigma^2 + tau^2))
    if(!unif) prior <- log(2) + dt_nonstandard(tau, sigma=25, log=TRUE) else prior <- log(1/50)
    return(exp(terms + prior))
}

tau <-  seq(0,50,len=200)
# note that these are not normalized and don't integrate to same value
plot(tau, sapply(tau, dtau, y, sigma, unif = FALSE), type = 'l')
lines(tau, sapply(tau, dtau, y, sigma, unif = TRUE), lty = 2)
legend('topright', legend = c('half-t', 'uniform'), lty = 1:2)

# normalized
d1 <- sapply(tau, dtau, y, sigma, unif = FALSE)
plot(tau, d1/sum(d1), type = 'l')
d2 <- sapply(tau, dtau, y, sigma, unif = TRUE)
lines(tau, d2/sum(d2), lty = 2)
legend('topright', legend = c('half-t', 'uniform'), lty = 1:2)


## @knitr jags-comparison

## jags comparison

library(R2jags)

# uniform prior
codeGelman <- function() {
    for(j in 1:J) {
        y[j] ~ dnorm(theta[j], isigma2[j])
        theta[j] ~ dnorm(mu, itau2)
    }
    mu ~ dnorm(0, .000001)
    itau2 <- 1/(tau*tau)
    tau ~ dunif(0,50)
}
out <- jags(data = list(y = y, J=J, isigma2 = 1/(sigma^2)), inits = list(list(tau = 10, theta = y, mu = mean(y))),
  parameters.to.save = c('theta','mu','tau'), n.chains = 1,
  n.iter = 15000, n.burnin = 1000, n.thin = 1, model.file = codeGelman, DIC = FALSE)
out.mcmc <- as.mcmc(out)[[1]]

ts.plot(out.mcmc[,'tau'])
hist(out.mcmc[,'tau'],breaks=30)

## @knitr schools-randomDist

## assess sensitivity to normality assumption on the random effects

codeGelmanT <- nimbleCode({
    for(j in 1:J) {
        y[j] ~ dnorm(theta[j], sd = sigma[j])
        theta[j] ~ dt(mu, sigma = tau, df = 5)
    }
    mu ~ dnorm(0, .000001)
    tau ~ dunif(0, 50)
})

modelGelmanT <- nimbleModel(codeGelmanT, data = list(y = y),
                                constants = list(sigma = sigma, J = J),
                                inits = list(theta = y, mu = mean(y),
                                             tau = sd(y)))

cmodelGelmanT <- compileNimble(modelGelmanT)
confGelmanT <- configureMCMC(modelGelmanT, monitors = c('theta', 'mu', 'tau'))
confGelmanT$removeSamplers('tau')
confGelmanT$addSampler('tau', type = 'RW', control = list(log = TRUE))
mcmcGelmanT <- buildMCMC(confGelmanT)
cmcmcGelmanT <- compileNimble(mcmcGelmanT, project = modelGelmanT)
cmcmcGelmanT$run(nIts)
smpGelmanT <- as.matrix(cmcmcGelmanT$mvSamples)[(burnin+1):nIts, ]

par(mfrow = c(1,2))
plot(density(smpGelman[[1]][ , 'mu']), xlim = c(-10, 30), xlab = expression(mu), main = '')
lines(density(smpGelmanT[ , 'mu']), lty = 2)
legend('topright', legend = c('normal', 't'), lty = 1:2)
plot(density(smpGelman[[1]][ , 'theta[1]']), xlim = c(-10, 50), xlab = expression(theta[1]), main = '')
lines(density(smpGelmanT[ , 'theta[1]']), lty = 2)
legend('topright', legend = c('normal', 't'), lty = 1:2)

par(mfrow = c(3,3))
itSample <- sample(1:(nIts-burnin), 9)
for(i in 1:9)
    hist(smpGelman[[1]][itSample[i] , 3:10], main = '', xlab = expression(theta))

## @knitr schools-postpred-check

## posterior-predictive checks of random effects distribution

# somewhat inefficient to do this as part of the MCMC vs. post-processing

codeGelmanPP <- nimbleCode({
    for(j in 1:J) {
        y[j] ~ dnorm(theta[j], sd = sigma[j])
        theta[j] ~ dnorm(mu, sd = tau)
        yrep[j] ~ dnorm(theta[j], sd = sigma[j])
    }
    mu ~ dnorm(0, .000001)
    tau ~ dunif(0, 50)
    
})

nIts <- 50000
modelGelmanPP <- nimbleModel(codeGelmanPP, data = list(y = y),
                                constants = list(sigma = sigma, J = J),
                                inits = list(theta = y, mu = mean(y),
                                             tau = sd(y)))

cmodelGelmanPP <- compileNimble(modelGelmanPP)
confGelmanPP <- configureMCMC(modelGelmanPP, monitors = c('yrep'))
confGelmanPP$removeSamplers('tau')
confGelmanPP$addSampler('tau', type = 'RW', control = list(log = TRUE))
mcmcGelmanPP <- buildMCMC(confGelmanPP)
cmcmcGelmanPP <- compileNimble(mcmcGelmanPP, project = modelGelmanPP)
cmcmcGelmanPP$run(nIts)
smpGelmanPP <- as.matrix(cmcmcGelmanPP$mvSamples)[(burnin+1):nIts, ]

skew <- function(x) {
    xbar <- mean(x)
    return(mean((x - xbar)^3) / (( mean((x - xbar)^2) )^(3/2)))
}

yrepMax <- apply(smpGelmanPP, 1, max)
yrepSkew <- apply(smpGelmanPP, 1, skew)

par(mfrow = c(1, 2))
hist(yrepMax, main = 'max yrep', xlab = 'max_yrep')
abline(v = max(y))
hist(yrepSkew, main = 'yrep skew', xlab = 'yrep skew')
abline(v = skew(y))

## @knitr schools-bayesFactor

joint_ytau <- function(tau, y, sigma, mu0, omega2, b) {
    J <- length(y)
    tau2 <- tau^2

    joint_y_onetau<- function(tau2, y, sigma, mu0, omega2, b) {
        v2js <- sigma^2 + tau2
        Vinv <- sum(1/v2js) + 1/omega2
        VinvM <- sum(y/v2js)+mu0/omega2
        M = VinvM/Vinv
        return((1/b)*exp(-(J/2)*log(2*pi) -0.5*sum(log(v2js)) -0.5*log(omega2)-0.5*log(Vinv) -
                    0.5*(sum(y^2/v2js)+mu0^2/omega2 - M*VinvM)))
    }

    out <- sapply(tau2, joint_y_onetau, y, sigma, mu0, omega2, b)
    return(out)
}

# check to see if unnormalized p(tau2|y) looks right
tau <- seq(0, 40, length = 100)
plot(tau, joint_ytau(tau, y, sigma, 0, 10^2, 50), type = 'l')


b <- 50; mu0 <- 0; omega2 <- 1^2
py1 <- integrate(joint_ytau, 0, b, y, sigma, mu0, omega2, b)

omega2 <- 10^2
py10 <- integrate(joint_ytau, 0, b, y, sigma, mu0, omega2, b)

omega2 <- 100^2
py100 <- integrate(joint_ytau, 0, b, y, sigma, mu0, omega2, b)

omega2 <- 10000^2
py10000 <- integrate(joint_ytau, 0, b, y, sigma, mu0, omega2, b)

b <- 500
py100b500 <- integrate(joint_ytau, 0, b, y, sigma, mu0, omega2, b)

# could actually consider different priors as different models
# and use the BF to compare them

# BF for 1 vs 10
py1$value / py10$value

# BF for 1 vs 100
py1$value / py100$value


## @knitr schools-bayesFactor-computation

# pretend integrals can't be done in closed-form
# and use the approach of sampling from the prior

K <- 100000
mu0 = 0
omega2 <- 10^2
b <- 50
tau <- runif(K, 0, b)
mu <- rnorm(K, mu0, sqrt(omega2))
theta <- matrix(rnorm(K*J, mu, tau), nrow = K)

yMat <- matrix(y, nrow = K, ncol = J, byrow = TRUE)
sigmaMat <- matrix(sigma, nrow = K, ncol = J, byrow = TRUE)
# careful here - exponentiation can lead to numerical underflow for larger datasets
logLik <- exp(rowSums(dnorm(yMat, theta, sigmaMat, log = TRUE)))
py10_sim <- mean(logLik)


omega2 <- 100^2
tau <- runif(K, 0, b)
mu <- rnorm(K, mu0, sqrt(omega2))
theta <- matrix(rnorm(K*J, mu, tau), nrow = K)

logLik <- exp(rowSums(dnorm(yMat, theta, sigmaMat, log = TRUE)))
py100_sim <- mean(logLik)

## @knitr schools-dic

schools_normal <- nimbleCode({
    for(j in 1:J) {
        y[j] ~ dnorm(theta[j], sd = sigma[j])
        theta[j] ~ dnorm(mu, sd = tau)
    }
    mu ~ dnorm(0, .000001)
    tau ~ dunif(0, 50)
})

schools_t <- nimbleCode({
    for(j in 1:J) {
        y[j] ~ dnorm(theta[j], sd = sigma[j])
        theta[j] ~ dt(mu, sigma = tau, df = 5)
    }
    mu ~ dnorm(0, .000001)
    tau ~ dunif(0, 50)
})


model_normal <- nimbleModel(schools_normal, data = list(y = y),
                                constants = list(sigma = sigma, J = J),
                                inits = list(theta = y, mu = mean(y),
                                             tau = sd(y)))
model_t <- nimbleModel(schools_t, data = list(y = y),
                                constants = list(sigma = sigma, J = J),
                                inits = list(theta = y, mu = mean(y),
                                             tau = sd(y)))
nIts <- 30000
burnin <- 1000

cmodel_normal <- compileNimble(model_normal)
conf_normal <- configureMCMC(model_normal, monitors = c('theta','mu','tau'))
conf_normal$removeSamplers('tau')
conf_normal$addSampler('tau', type = 'RW', control = list(log = TRUE))
mcmc_normal <- buildMCMC(conf_normal)
cmcmc_normal <- compileNimble(mcmc_normal, project = model_normal)
cmcmc_normal$run(nIts)
smp_normal <- as.matrix(cmcmc_normal$mvSamples)[(burnin+1):nIts, ]

cmodel_t <- compileNimble(model_t)
conf_t <- configureMCMC(model_t, monitors = c('theta','mu','tau'))
conf_t$removeSamplers('tau')
conf_t$addSampler('tau', type = 'RW', control = list(log = TRUE))
mcmc_t <- buildMCMC(conf_t)
cmcmc_t <- compileNimble(mcmc_t, project = model_t)
cmcmc_t$run(nIts)
smp_t <- as.matrix(cmcmc_t$mvSamples)[(burnin+1):nIts, ]

dicFull <- function(y, sigma, thetaSmp) {
    yMat <- matrix(y, nrow = nrow(thetaSmp), ncol = length(y), byrow = TRUE)
    sigmaMat <- matrix(sigma, nrow = nrow(thetaSmp), ncol = length(sigma), byrow = TRUE)
    deviance <- rowSums(-2*dnorm(yMat, thetaSmp, sigmaMat, log = TRUE))
    dbar <- mean(deviance)
    dThetaBar <- -2*sum(dnorm(y, colMeans(thetaSmp), sigma, log = TRUE))
    return(c(dic = 2*dbar - dThetaBar, p = dbar-dThetaBar))
}

dicMarginal <- function(y, sigma, phiSmp) {
    J <- length(y)
    local <- function(phi) {
        covar <- matrix(phi['tau']^2, J, J)
        diag(covar) <- diag(covar) + sigma^2
        ch <- chol(covar)
        return(dmnorm_chol(y, phi['mu'], ch, prec_param = FALSE, log = TRUE))
    }
    deviance <- -2*apply(phiSmp, 1, local)
    dbar <- mean(deviance)
    dThetaBar <- -2*local(colMeans(phiSmp))
    return(c(dic = 2*dbar - dThetaBar, p = dbar-dThetaBar))
}

# DIC values for normal model
dicFull(y, sigma, smp_normal[ , 3:(J+2)])
dicMarginal(y, sigma, smp_normal[ , 1:2])
        
# DIC values for t model
dicFull(y, sigma, smp_t[ , 3:(J+2)])
dicMarginal(y, sigma, smp_t[ , 1:2])
       
## @knitr schools-waic

# could write this algorithm as a nimbleFunction for greater efficiency and generality

waic <- function(y, sigma, thetaSmp) {
    y <- matrix(y, nrow = nrow(thetaSmp), ncol = length(y), byrow = TRUE)
    sigma <- matrix(sigma, nrow = nrow(thetaSmp), ncol = length(sigma), byrow = TRUE)
    lik <- dnorm(y, thetaSmp, sigma)
    lppd <- sum(log(colMeans(lik)))
    pwaic2 <- sum(apply(log(lik), 2, var))
    return(c(waic = -2*lppd + 2*pwaic2, p = pwaic2))
}

# WAIC for normal model
waic(y, sigma, smp_normal)
# WAIC for t model
waic(y, sigma, smp_t)

