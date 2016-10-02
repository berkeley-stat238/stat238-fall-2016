## @knitr bliss-bclt

w <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
y <- c(6, 13, 18, 28, 52, 53, 61, 60)
n <- c(59, 60, 62, 56, 63, 59, 62, 60)

logLik <- function(theta){
  m1 <- exp(theta[3])
  sigma <- exp(theta[2])
  mu <- theta[1]
  x <- (w-mu)/sigma
  logp <- m1*(x-log(1+exp(x)))
  return(sum(y*logp+(n-y)*log(1-exp(logp))))
}

dinvgamma <- function(x, alpha, beta, log = TRUE) {
    dens <- alpha*log(beta) - (alpha + 1)*log(x) -beta/x -lgamma(alpha)
    if(log) return(dens) else return(exp(dens))
}
    
logPrior <- function(theta, hypers){
  m1 <- exp(theta[3])
  sigma2 <- (exp(theta[2]))^2
  mu <- theta[1]
  return(dgamma(m1, hypers$a0, hypers$b0, log = TRUE) +
         dnorm(mu, hypers$c0, hypers$d0, log = TRUE) +
         dinvgamma(sigma2, hypers$e0, hypers$f0))
}

logPriorTrans <- function(theta, hypers) {
    mu <- theta[1]
    logm1 <- theta[3]
    logsigma <- theta[2]
    sigma2 <- exp(2*logsigma)
    m1 <- exp(logm1)
    # log(2) irrelevant but here for clarity
    return(dnorm(mu, hypers$c0, hypers$d0, log = TRUE) +
           # include jacobian of inverse transformation in each of next two lines
           dinvgamma(sigma2, hypers$e0, hypers$f0, log = TRUE) + log(2) + 2*logsigma +
           dgamma(m1, hypers$a0, hypers$b0, log = TRUE) + logm1)
}
   

hypers <- list(a0 = .25, b0 = .25, c0 = 2, d0 = 10, e0 = 2, f0 = .001)

# get MLE

# starting values for optimization
theta <- c(1.75, log(0.1), log(1))

expit <- function(x) exp(x)/(1+exp(x))
x <- (w - theta[1]) / exp(theta[2])
p <- (expit(x))^exp(theta[3])
plot(w, y/n)
points(w, p, col = 'red')

out <- optim(theta, logLik, control = list(fnscale = -1))
theta <- out$par

post <- function(theta) logLik(theta) + logPriorTrans(theta, hypers)
out <- optim(theta, post, control = list(fnscale = -1), hessian = TRUE)

# BCLT

post_mode <- out$par
post_cov <- -solve(out$hessian)
post_sds <- sqrt(diag(post_cov))

p <- length(post_mode)
m <- 500

bclt_samples <- t(post_mode + t(chol(post_cov)) %*% matrix(rnorm(m*p), nrow = p))

# functional of interest
dose <- 1.8
pdose <- function(theta, dose) {
  m1 <- exp(theta[3])
  sigma <- exp(theta[2])
  mu <- theta[1]
  x <- (dose-mu)/sigma
  logp <- m1*(x-log(1+exp(x)))
  return(exp(logp))
}

pred_bclt <- apply(bclt_samples, 1, pdose, dose)
mean(pred_bclt)
sd(pred_bclt)


## @knitr bliss-grid

# grid-based approximation

lowers <- post_mode - 5*post_sds
uppers <- post_mode + 5*post_sds

k <- 100 # discretization

mu_grd <- seq(lowers[1], uppers[1], length = k)
logsigma_grd <- seq(lowers[2], uppers[2], length = k)
logm1_grd <- seq(lowers[3], uppers[3], length = k)
grd <- expand.grid(mu_grd, logsigma_grd, logm1_grd)


postVals <- apply(grd, 1, post)
postVals_prob <- exp(postVals) / sum(exp(postVals))

ids <- sample(seq_along(postVals_prob), size = m, replace = TRUE,
              prob = postVals_prob)
grid_samples <- grd[ids, ]

pred_grid <- apply(grid_samples, 1, pdose, dose)
mean(pred_grid)
sd(pred_grid)

plot(exp(bclt_samples[,2]), exp(bclt_samples[,3]))
plot(exp(grid_samples[,2]), exp(grid_samples[,3]))

plot(bclt_samples[,2], bclt_samples[,3])
plot(grid_samples[,2], grid_samples[,3])


## @knitr bliss-is

library(mvtnorm)
library(methods)  # needed when compiling within Sweave

Sis <- 1000000
gdraws_uncentered <- rmvt(Sis, sigma = post_cov, df = 1)
gdraws <- t(post_mode + t(gdraws_uncentered))

qLogDensity <- apply(gdraws, 1, post)
qLogDensity[is.na(qLogDensity)] <- -Inf
gLogDensity <- dmvt(gdraws_uncentered, sigma = post_cov, df = 1, log = TRUE)
wgt <- exp(qLogDensity - gLogDensity)    
wgtMn <- mean(wgt, na.rm = TRUE)

# Gelman (10.4)
seff <- 1/sum( (wgt/sum(wgt)) ^2)  

# use IS for response at dose = 1.8
pm <- mean(apply(gdraws, 1, pdose, dose)*wgt, na.rm =TRUE) / wgtMn
p2mom <- mean(apply(gdraws, 1, pdose, dose)^2*wgt, na.rm = TRUE) / wgtMn
psd <- sqrt(p2mom - pm^2)


# samples based on IS approach - SIR

ids <- sample(seq_len(Sis), m, replace = TRUE, prob = wgt)
is_samples <- gdraws[ids, ]

plot(exp(is_samples[,2]), exp(is_samples[,3]))
plot(is_samples[,2], is_samples[,3])

## @knitr schools-gibbs

set.seed(1)
J <- 8
n <- rep(40, J)
sigma <- 5
tau <- 5  # 0.5 is case where mixing is not good
mu <- 1
y <- matrix(0, nrow = J, ncol = max(n))
theta <- rnorm(J, mu, tau)

for(j in 1:J)
    y[j, 1:n[j]] <- rnorm(n[j], theta[j], sigma)

code <- nimbleCode({
    itau2 ~ dgamma(a_tau, b_tau)
    isigma2 ~ dgamma(a_sigma, b_sigma)
    for(j in 1:J) {
        for(i in 1:n[j])
            y[j, i] ~ dnorm(theta[j], isigma2)
        theta[j] ~ dnorm(mu, itau2)
    }
    mu ~ dnorm(0, .000001)
})

eps = .001
m <- nimbleModel(code, data = list(y = y),
                 constants = list(n = n, J = J, a_tau = eps,
                                  b_tau = eps, a_sigma = eps, b_sigma =eps),
                 inits = list(itau2 = 1/8, isigma2 =1/8,
                              theta = rep(0, J), mu = 0))


conf <- configureMCMC(m, monitors = c('itau2', 'isigma2', 'mu',
                                      'theta'))

cm <- compileNimble(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(10000)

smp <- as.matrix(cmcmc$mvSamples)

sigma <- sqrt(1/smp[ , 'isigma2'])
tau <- sqrt(1/smp[ , 'itau2'])
theta1 <- smp[ , 'theta[1]']

par(mfrow = c(2,3))
ts.plot(sigma)
ts.plot(tau)
ts.plot(smp[ , 'mu'])
plot(sigma, tau)
ts.plot(theta1)
plot(tau, theta1)

cor(tau[100:1000], theta1[100:1000])

## @knitr bliss-mh


dloggamma <- nimbleFunction(
    run = function(x = double(0), a = double(0), b = double(0),
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- a*x-b*exp(x)
        if(log) return(logProb)
        else return(exp(logProb)) 
        })

# WARNING: just a placeholder we know is not used in NIMBLE's Metropolis-Hastings sampler
rloggamma <- nimbleFunction(
    run = function(n = integer(0), a = double(0), b = double(0)) {
        returnType(double(0))
        if(n != 1) print("rloggamma only allows n = 1; using n = 1.")
        return(0)
    })


dlogsqrtIG <- nimbleFunction(
    run = function(x = double(0), a = double(0), b = double(0),
        log = integer(0, default = 0)) {
        returnType(double(0))
        logProb <- -2*a*x - b*exp(-2*x)
        if(log) return(logProb)
        else return(exp(logProb)) 
        })

# WARNING: just a placeholder we know is not used in NIMBLE's Metropolis-Hastings sampler
rlogsqrtIG <- nimbleFunction(
    run = function(n = integer(0), a = double(0), b = double(0)) {
        returnType(double(0))
        if(n != 1) print("rlogsqrtIG only allows n = 1; using n = 1.")
        return(0)
    })

registerDistributions(list(
        dloggamma = list(
               BUGSdist = "dloggamma(a,b)"),
        dlogsqrtIG = list(
               BUGSdist = "dlogsqrtIG(a,b)")))


code <- nimbleCode({
    for(i in 1:G) {
        x[i] <- (w[i] - mu) / sigma
        p[i] <- (expit(x[i]))^m1
        y[i] ~ dbin(p[i], n[i])
    }
    m1 <- exp(logm1)
    logm1 ~ dloggamma(a0, b0)
    mu ~ dnorm(c0, sd = d0)
    logsigma ~ dlogsqrtIG(e0, f0)
    sigma <- exp(logsigma)
})

inits <- list(mu = 2, logsigma = -8, logm1 = -4)
m <- nimbleModel(code, data = list(w = w, y = y),
                 constants = c(list(n = n, G = length(n)), hypers),
                 inits = inits)

params <- c('mu','logsigma','logm1')

conf <- configureMCMC(m)
conf$getSamplers()
conf$removeSamplers(params)
conf$addSampler(params, type = 'RW_block',
                control = list(adaptive = FALSE, scale = 1,
                               propCov = diag(c(0.00012, 0.33, .10))))
# the proposal variances are roughly the posterior variances

cm <- compileNimble(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
nIts <- 4000
set.seed(0)
cmcmc$run(nIts)

smp <- as.matrix(cmcmc$mvSamples)

par(mfrow = c(2,3))
ts.plot(smp[ , 'mu'])
ts.plot(smp[ , 'logsigma'])
ts.plot(smp[ , 'logm1'])
plot(smp[ , 'mu'], smp[ , 'logsigma'])
plot(smp[ , 'mu'], smp[ , 'logm1'])
plot(smp[ , 'logsigma'], smp[ , 'logm1'])

# save for later for convergence assessment
cmA <- cm
cmcmcA <- cmcmc
smpA1 <- smp[ , params] 


## @knitr bliss-mh-corr

# ordering of params in output not the same as when use addSampler - careful

c2 <- (2.4/sqrt(length(params)))^2 # theoretical result Gelman p. 296
propCov <- c2 * cov(smp[500:nIts, params])

conf <- configureMCMC(m)
conf$getSamplers()
conf$removeSamplers(params)
conf$addSampler(params, type = 'RW_block',
                control = list(adaptive = FALSE, scale = 1,
                               propCov = propCov))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)

inits <- list(mu = 2, logsigma = -8, logm1 = -4)
cm$setInits(inits)

set.seed(0)
cmcmc$run(nIts)

smp <- as.matrix(cmcmc$mvSamples)

par(mfrow = c(2,3))
ts.plot(smp[ , 'mu'])
ts.plot(smp[ , 'logsigma'])
ts.plot(smp[ , 'logm1'])
plot(smp[ , 'mu'], smp[ , 'logsigma'])
plot(smp[ , 'mu'], smp[ , 'logm1'])
plot(smp[ , 'logsigma'], smp[ , 'logm1'])

# acceptance rate:
mean(smp[2:nIts,1]!=smp[1:(nIts-1),1])

# careful of ordering of params in output matrix
smp_pdose <- apply(smp[ , params], 1, pdose, dose = dose)
mean(smp_pdose)
sd(smp_pdose)

# save for later for convergence assessment
cmB <- cm
cmcmcB <- cmcmc
smpB1 <- smp[ , params]

## @knitr bliss-conv

# let's re-run our two samplers with multiple starting points

# different starting values 

inits <- list(mu = 1.5, logsigma = -2, logm1 = 1)
cmA$setInits(inits)

set.seed(0)
cmcmcA$run(nIts)
smpA2 <- as.matrix(cmcmcA$mvSamples)[ , params]

cmB$setInits(inits)

set.seed(0)
cmcmcB$run(nIts)
smpB2 <- as.matrix(cmcmcB$mvSamples)[ , params]

# traceplots

par(mfrow = c(2,3))
ts.plot(smpA1[ , 'mu'])
ts.plot(smpA1[ , 'logsigma'])
ts.plot(smpA1[ , 'logm1'])
ts.plot(smpB[ , 'mu'])
ts.plot(smpB[ , 'logsigma'])
ts.plot(smpB[ , 'logm1'])

# effective sample size
nonBurn <- 2001:nIts   # note not needed for 2nd set as it started 'warmed up'
# but apply to both for comparability

library(coda)
apply(smpA1[nonBurn, ], 2, effectiveSize)
apply(smpB[nonBurn, ], 2, effectiveSize)

# ACF
apply(smpA1[nonBurn, ], 2, acf)
apply(smpB[nonBurn, ], 2, acf)

# potential scale reduction factor (Gelman-Rubin)
firstHalf <- 2001:3000
secondHalf <- 3001:4000
mcmcList <- mcmc.list(as.mcmc(smpA1[firstHalf, ]), as.mcmc(smpA1[secondHalf, ]),
                      as.mcmc(smpA2[firstHalf, ]), as.mcmc(smpA2[secondHalf, ]))
plot(mcmcList)

                      
gelman.diag(mcmcList)  # not good!

mcmcList <- mcmc.list(as.mcmc(smpB1[firstHalf, ]), as.mcmc(smpB1[secondHalf, ]),
                      as.mcmc(smpB2[firstHalf, ]), as.mcmc(smpB2[secondHalf, ]))
gelman.diag(mcmcList)  # better

## @knitr bliss-mh-adaptive

### univariate

conf <- configureMCMC(m)
conf$getSamplers()

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
# allow time for adaptation
nIts <- 5000
cmcmc$run(nIts)

smp <- as.matrix(cmcmc$mvSamples)

par(mfrow = c(2,3))
ts.plot(smp[ , 'mu'])
ts.plot(smp[ , 'logsigma'])
ts.plot(smp[ , 'logm1'])
plot(smp[ , 'mu'], smp[ , 'logsigma'])
plot(smp[ , 'mu'], smp[ , 'logm1'])
plot(smp[ , 'logsigma'], smp[ , 'logm1'])

# acceptance rate:
mean(smp[2001:nIts,1]!=smp[2000:(nIts-1),1])


### multivariate

conf <- configureMCMC(m)
conf$getSamplers()
conf$removeSamplers(params)
# need to start at reasonable univariate scales or learning of
# covariance can take a long time because of terrible initial mixing
conf$addSampler(params, type = 'RW_block',
                control = list(propCov = diag(c(0.00012, 0.33, .10))))

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m, resetFunctions = TRUE)
# allow time for adaptation
nIts <- 5000
cmcmc$run(nIts)

smp <- as.matrix(cmcmc$mvSamples)

par(mfrow = c(2,3))
ts.plot(smp[ , 'mu'])
ts.plot(smp[ , 'logsigma'])
ts.plot(smp[ , 'logm1'])
plot(smp[ , 'mu'], smp[ , 'logsigma'])
plot(smp[ , 'mu'], smp[ , 'logm1'])
plot(smp[ , 'logsigma'], smp[ , 'logm1'])

# acceptance rate:
mean(smp[2001:nIts,1]!=smp[2000:(nIts-1),1])


## @knitr schools-reparam

set.seed(1)
J <- 8
n <- rep(40, J)
sigma <- 5
tau <- 0.5  # 0.5 is case where mixing is not good
mu <- 1
y <- matrix(0, nrow = J, ncol = max(n))
theta <- rnorm(J, mu, tau)

for(j in 1:J)
    y[j, 1:n[j]] <- rnorm(n[j], theta[j], sigma)

code <- nimbleCode({
    logalpha ~ dnorm(0, .0001)
    tau <- sqrt(1/itau2)
    itau2 ~ dgamma(a_tau, b_tau)
    isigma2 ~ dgamma(a_sigma, b_sigma)
    for(j in 1:J) {
        for(i in 1:n[j])
            y[j, i] ~ dnorm(mu + tau*theta[j], isigma2)
        theta[j] ~ dnorm(0, 1)
    }
    mu ~ dnorm(0, .000001)
})

eps = .001
m <- nimbleModel(code, data = list(y = y),
                 constants = list(n = n, J = J, a_tau = eps,
                                  b_tau = eps, a_sigma = eps, b_sigma =eps),
                 inits = list(itau2 = 1/8, isigma2 =1/8,
                              theta = rep(0, J), mu = 0))


conf <- configureMCMC(m, monitors = c('itau2', 'isigma2', 'mu',
                                      'theta'))
conf$removeSamplers('itau2')
conf$addSampler('itau2', 'RW', control = list(log = TRUE))

cm <- compileNimble(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(10000)

smp <- as.matrix(cmcmc$mvSamples)

sigma <- sqrt(1/smp[ , 'isigma2'])
tau <- sqrt(1/smp[ , 'itau2'])
theta1 <- smp[ , 'theta[1]']
mu <- smp[ , 'mu']
par(mfrow = c(2,3))
ts.plot(sigma)
ts.plot(tau)
ts.plot(mu)
plot(sigma, tau)
ts.plot(mu + tau*theta1)
plot(tau, theta1)

## @knitr bliss-laplace

logLik <- function(theta){
  m1 <- exp(theta[3])
  sigma <- exp(theta[2])
  mu <- theta[1]
  x <- (w-mu)/sigma
  logp <- m1*(x-log(1+exp(x)))
  return(sum(y*logp+(n-y)*log(1-exp(logp))))
}

dinvgamma <- function(x, alpha, beta, log = TRUE) {
    dens <- alpha*log(beta) - (alpha + 1)*log(x) -beta/x -lgamma(alpha)
    if(log) return(dens) else return(exp(dens))
}
    
logPriorTrans <- function(theta, hypers) {
    mu <- theta[1]
    logm1 <- theta[3]
    logsigma <- theta[2]
    sigma2 <- exp(2*logsigma)
    m1 <- exp(logm1)
    # log(2) irrelevant but here for clarity
    return(dnorm(mu, hypers$c0, hypers$d0, log = TRUE) +
           # include jacobian of inverse transformation in each of next two lines
           dinvgamma(sigma2, hypers$e0, hypers$f0, log = TRUE) + log(2) + 2*logsigma +
           dgamma(m1, hypers$a0, hypers$b0, log = TRUE) + logm1)
}

dose <- 1.8
pdose <- function(theta, dose) {
  m1 <- exp(theta[3])
  sigma <- exp(theta[2])
  mu <- theta[1]
  x <- (dose-mu)/sigma
  logp <- m1*(x-log(1+exp(x)))
  return(exp(logp))
}

hypers <- list(a0 = .25, b0 = .25, c0 = 2, d0 = 10, e0 = 2, f0 = .001)
theta <- c(1.75, log(0.1), log(1))


post <- function(theta) logLik(theta) + logPriorTrans(theta, hypers)
L <- optim(theta, post, control = list(fnscale = -1), hessian = TRUE)

# posterior mean

h <- function(theta, dose) logLik(theta) + logPriorTrans(theta, hypers) + log(pdose(theta, dose))
Lstar <- optim(theta, h, control = list(fnscale = -1), hessian = TRUE, dose = dose)

SigmaStar <- -solve(Lstar$hessian)
Sigma <- -solve(L$hessian)
logSqrtDetermStar <- sum(log(diag(chol(SigmaStar))))
logSqrtDeterm <- sum(log(diag(chol(Sigma))))
pm <- exp(logSqrtDetermStar - logSqrtDeterm + Lstar$value - L$value)

# finding the explicit inverse is not the recommended numerical method; for larger matrices, make use of Cholesky decomposition and use backsolves
        

# posterior 2nd moment

h <- function(theta, dose) logLik(theta) + logPriorTrans(theta, hypers) + 2*log(pdose(theta, dose))
Lstar <- optim(theta, h, control = list(fnscale = -1), hessian = TRUE, dose = dose)

SigmaStar <- -solve(Lstar$hessian)
logSqrtDetermStar <- sum(log(diag(chol(SigmaStar))))

p2mom <- exp(logSqrtDetermStar - logSqrtDeterm + Lstar$value - L$value)
psd <- sqrt(p2mom - pm^2)

# comparison with BCLT

pm <- pdose(L$par, dose = dose)


## @knitr misc

#### Below not used

## @knitr schools-t

set.seed(1)
J <- 8
n <- rep(40, J)
sigma <- 5
tau <- 5  # 0.5 is case where mixing is not good
mu <- 1
y <- matrix(0, nrow = J, ncol = max(n))
theta <- rnorm(J, mu, tau)

for(j in 1:J)
    y[j, 1:n[j]] <- rnorm(n[j], theta[j], sigma)

code <- nimbleCode({
    itau2 ~ dgamma(a_tau, b_tau)
    isigma2 ~ dgamma(a_sigma, b_sigma)
    for(j in 1:J) {
        for(i in 1:n[j])
            y[j, i] ~ dt(theta[j], isigma2, df = nu)
        theta[j] ~ dnorm(mu, itau2)
    }
    mu ~ dnorm(0, .000001)
})

eps = .001
m <- nimbleModel(code, data = list(y = y),
                 constants = list(n = n, J = J, a_tau = eps,
                                  b_tau = eps, a_sigma = eps, b_sigma =eps,
                                  nu = 3),
                 inits = list(itau2 = 1/8, isigma2 =1/8,
                              theta = rep(0, J), mu = 0))


conf <- configureMCMC(m, monitors = c('itau2', 'isigma2', 'mu',
                                      'theta'))
conf$printSamplers()

cm <- compileNimble(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(1000)

smp <- as.matrix(cmcmc$mvSamples)

sigma <- sqrt(1/smp[ , 'isigma2'])
tau <- sqrt(1/smp[ , 'itau2'])
theta1 <- smp[ , 'theta[1]']

par(mfrow = c(2,3))
ts.plot(sigma)
ts.plot(tau)
ts.plot(smp[ , 'mu'])
plot(sigma, tau)
ts.plot(theta1)
plot(tau, theta1)


## @knitr schools-cent-uncent

set.seed(1)
J <- 8
n <- rep(40, J)
sigma <- 5
tau <- 5
mu <- 1
y <- matrix(0, nrow = J, ncol = max(n))
theta <- rnorm(J, mu, tau)

for(j in 1:J)
    y[j, 1:n[j]] <- rnorm(n[j], theta[j], sigma)

code <- nimbleCode({
    itau2 ~ dgamma(a_tau, b_tau)
    isigma2 ~ dgamma(a_sigma, b_sigma)
    for(j in 1:J) {
        for(i in 1:n[j])
            y[j, i] ~ dnorm(mu + theta[j], isigma2)
        theta[j] ~ dnorm(0, itau2)
    }
    mu ~ dnorm(0, .000001)
})

eps = .001
m <- nimbleModel(code, data = list(y = y),
                 constants = list(n = n, J = J, a_tau = eps,
                                  b_tau = eps, a_sigma = eps, b_sigma =eps),
                 inits = list(itau2 = 1/8, isigma2 =1/8,
                              theta = rep(0, J), mu = 0))


conf <- configureMCMC(m, monitors = c('itau2', 'isigma2', 'mu',
                                      'theta'))

cm <- compileNimble(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(10000)

smp <- as.matrix(cmcmc$mvSamples)

sigma <- sqrt(1/smp[ , 'isigma2'])
tau <- sqrt(1/smp[ , 'itau2'])
theta1 <- smp[ , 'theta[1]']

par(mfrow = c(2,3))
ts.plot(sigma)
ts.plot(tau)
ts.plot(smp[ , 'mu'])
plot(sigma, tau)
ts.plot(theta1)
plot(tau, theta1)


## @knitr schools-t-augment

# note that to keep conjugacy 
code <- nimbleCode({
    itau2 ~ dgamma(a_tau, b_tau)
    isigma2 ~ dgamma(a_sigma, b_sigma)
    for(j in 1:J) {
        for(i in 1:n[j]) {
            y[j, i] ~ dnorm(theta[j], q[j, i])
            q[j, i] ~ dgamma(nu/2, nu / (2*isigma2))
        }
        theta[j] ~ dnorm(mu, itau2)
    }
    mu ~ dnorm(0, .000001)
})

eps = .001
m <- nimbleModel(code, data = list(y = y),
                 constants = list(n = n, J = J, a_tau = eps,
                                  b_tau = eps, a_sigma = eps, b_sigma =eps,
                                  nu = 3),
                 inits = list(itau2 = 1/8, isigma2 =1/8,
                              theta = rep(0, J), mu = 0))


conf <- configureMCMC(m, monitors = c('itau2', 'isigma2', 'mu',
                                      'theta'))
conf$printSamplers()

cm <- compileNimble(m)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)
cmcmc$run(1000)

smp <- as.matrix(cmcmc$mvSamples)

sigma <- sqrt(1/smp[ , 'isigma2'])
tau <- sqrt(1/smp[ , 'itau2'])
theta1 <- smp[ , 'theta[1]']

par(mfrow = c(2,3))
ts.plot(sigma)
ts.plot(tau)
ts.plot(smp[ , 'mu'])
plot(sigma, tau)
ts.plot(theta1)
plot(tau, theta1)

## @knitrs schools-vb

y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)
K <- length(y)

set.seed(0)
nIts <- 100

Malpha <- S2alpha <- matrix(0, nrow = nIts, ncol = K)
Mmu <- S2mu <- btau2 <- rep(0, nIts)
Malpha[1, ] <- rnorm(K)
S2alpha[1, ] <- runif(K)
Mmu[1] <- rnorm(1)
S2mu[1] <- runif(1)
btau2[1] <- runif(1)
atau2 <- (K-1)/2  # not a function of other variational parameters, so fixed

for(t in 2:nIts) {
    E1tau2 <- atau2/btau2[t-1]
    # update alpha variational parameters as a group
    Malpha[t, ] <- (y/sigma^2 + E1tau2*Mmu[t-1])  / (1/sigma^2 + E1tau2)
    S2alpha[t, ] <- 1 / (1/sigma^2 + E1tau2)
    # update mu variational parameters
    Mmu[t] <- mean(Malpha[t, ])
    S2mu[t] <- 1 / (K * E1tau2)
    # update tau2 variational parameter
    btau2[t] <- 0.5*(sum(S2alpha[t, ]) + sum(Malpha[t, ]^2) - 2*sum(Malpha[t, ])*Mmu[t] + K*Mmu[t]^2 + K*S2mu[t])
}

par(mfrow = c(1, 3))
ts.plot(Mmu)
ts.plot(S2mu)
ts.plot(sqrt((2/7)*btau2))  # 2/7 converts to BDA scaling and sqrt to BDA Fig 13.5

# variational approximation for tau

tau2 <- 1/ rgamma(10000, atau2, btau2[nIts])
par(mfrow = c(1, 2))
plot(density(sqrt(tau2)), xlim = c(0, 30),
     xlab = expression(tau), ylab = expression(g(tau,y)))

# note mismatch compared to exact result
min(sqrt(tau2))

# exact marginal posterior for tau2
# see BDA p. 117
tauPost <- function(tau, y, sigma) {
    tau2 <- tau^2
    vhatinv <- sum(1/(sigma^2 + tau2))
    muhat <- sum(y/(sigma^2 + tau2)) / vhatinv
    return(exp(-0.5*log(vhatinv) -0.5 * sum(log(sigma^2 + tau2)) - 0.5*sum((y - muhat)^2 / (sigma^2 + tau2))))
}

tau <- seq(0, 30, length = 100)

plot(tau, sapply(tau, tauPost, y, sigma), type = 'l',
     xlab = expression(tau), ylab = expression(p(tau, y)))

    
## @knitrs schools-vb-nonconj


galpha_fun <- function(x, y, sigma, Emu, E1tau2) {
    exp(-0.5*(x^2 * (1/sigma^2 + E1tau2) - 2*x*(y/sigma^2 +E1tau2*Emu)))
}
Ealpha_fun <- function(x, y, sigma, Emu, E1tau2, const) {
    const * x * exp(-0.5*(x^2 * (1/sigma^2 + E1tau2) - 2*x*(y/sigma^2 +E1tau2*Emu)))
}
Ealpha2_fun <- function(x, y, sigma, Emu, E1tau2, const) {
    const * x^2 * exp(-0.5*(x^2 * (1/sigma^2 + E1tau2) - 2*x*(y/sigma^2 +E1tau2*Emu)))
}

gmu_fun <- function(x, K, Ealpha, E1tau2) {
    exp(-0.5*E1tau2 * (sum(Ealpha^2) -2*x*sum(Ealpha) + K*x^2))
}
Emu_fun <- function(x, K, Ealpha, E1tau2, const) {
    const*x*exp(-0.5*E1tau2 * (sum(Ealpha^2) -2*x*sum(Ealpha) + K*x^2))
}
Emu2_fun <- function(x, K, Ealpha, E1tau2, const) {
    const*x^2*exp(-0.5*E1tau2 * (sum(Ealpha^2) -2*x*sum(Ealpha) + K*x^2))
}

gtau_fun <- function(x, K, Ealpha, Ealpha2, Emu, mu2) {
    exp(-K*log(x) -(0.5/x^2)*sum(Ealpha2 - 2*Ealpha*Emu + Emu2))
}
E1tau2_fun <- function(x, K, Ealpha, Ealpha2, Emu, mu2, const) {
    const*(1/x^2) * exp(-K*log(x) -(0.5/x^2)*sum(Ealpha2 - 2*Ealpha*Emu + Emu2))
}

set.seed(0)
nIts <- 100

Ealpha <- rnorm(K)
E1tau2 <- runif(1)
Emu <- rnorm(1)
# ensure constraints on 2nd moment given squared first moment
Ealpha2 <- Ealpha^2 + runif(K)
Emu2 <- Emu^2 + rnorm(1)
calpha <- rep(0, K)

for(t in 2:nIts) {
    for(k in 1:K) {
        calpha[k] <- integrate(galpha_fun, lower = -Inf, upper = Inf, y[k], sigma[k], Emu, E1tau2)$value
        Ealpha[k] <- integrate(Ealpha_fun, lower = -Inf, upper = Inf, y[k], sigma[k], Emu, E1tau2, 1/calpha[k])$value
        Ealpha2[k] <- integrate(Ealpha2_fun, lower = -Inf, upper = Inf, y[k], sigma[k], Emu, E1tau2, 1/calpha[k])$value
    }
    cmu <- integrate(gmu_fun, lower = -Inf, upper = Inf, K, Ealpha, E1tau2)$value
    Emu <- integrate(Emu_fun, lower = -Inf, upper = Inf, K, Ealpha, E1tau2, 1/cmu)$value
    Emu2 <- integrate(Emu2_fun, lower = -Inf, upper = Inf, K, Ealpha, E1tau2, 1/cmu)$value
    ctau <- integrate(gtau_fun, lower = 0, upper = Inf, K, Ealpha, Ealpha2, Emu, Emu2)$value
    E1tau2 <- integrate(E1tau2_fun, lower = 0, upper = Inf, K, Ealpha, Ealpha2, Emu, Emu2, 1/ctau)$value
}

# density for tau2

plot(tau, (1/ctau)*gtau_fun(tau, K, Ealpha, Ealpha2, Emu, Emu2), type = 'l',
     xlab = expression(tau), ylab = expression(p(tau, y)))

# credible interval (not clear we'd want this given approximation involved in VB...)
# manual 'optimization' to invert the CDF
(1/ctau) * integrate(gtau_fun, 0, 5.9, K, Ealpha, Ealpha2, Emu, Emu2)$value
(1/ctau) * integrate(gtau_fun, 18.5, Inf, K, Ealpha, Ealpha2, Emu, Emu2)$value
