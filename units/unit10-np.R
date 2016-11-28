## @knitr bias-variance

par(mfrow = c(2, 2), mgp = c(1.8, .7, 0), mai = c(.3, .3, .3, .3))

x <- seq(0, 1, len = 40)
f <- sin(2*pi*x)
y <- rnorm(40, f, .3)
plot(x, y, ylim = c(-1.5, 1.5))
lines(x, f, col = 2)
# fit to one simulated dataset
mod <- ksmooth(x, y, kernel = "normal", x.points = x, bandwidth = .25)
lines(mod$x, mod$y)
title('oversmoothing')
print(sum((mod$y-f)^2))
legend('topright', col = c('black','red'), legend = c('fit', 'truth'),
       lty = rep(1, 2))

plot(x, f, col = 2, type = 'l', ylim = c(-1.5, 1.5))
# fit to many simulated datasets
for(i in 1:100){
  y <- rnorm(40, f, .3)
  mod <- ksmooth(x, y, kernel = "normal", bandwidth = .25)
  lines(mod$x, mod$y, lwd = .5)
}
lines(x, f, col = 2)

plot(x, y, ylim = c(-1.5, 1.5))
lines(x, f, col = 2)
# fit to one simulated dataset
mod <- ksmooth(x, y, kernel = "normal", x.points = x, bandwidth = .05)
lines(mod$x, mod$y)
print(sum((mod$y-f)^2))
title('undersmoothing')

# fit to many simulated datasets
x <- seq(0, 1, len = 40)
f <- sin(2*pi*x)
plot(x, f, col = 2, type = 'l', ylim = c(-1.5, 1.5))
for(i in 1:100){
  y <- rnorm(40, f, .3)
  mod <- ksmooth(x, y, kernel = "normal", bandwidth = .05)
  lines(mod$x, mod$y, lwd = .5)
}
lines(x, f, col = 2)


## @knitr bspline

library(splines)

x <- seq(0, 1, len = 100)
kn <- seq(.1, .9, len = 6)
B <- bs(x, knots = kn)

par(mfrow = c(1, 1),mgp = c(1.8, .7, 0), mai = c(.3, .3, .3, .3))
plot(x, B[, 1], ylim = c(-1.5, 1.5), type = 'n')
for(i in 1:9){
  lines(x, B[, i], col = i)
}

set.seed(0)
betas <- rnorm(9)

f <- B%*%betas
lines(x, f, lwd = 2)

betas

vals <- c(.07, .15, .25, .4, .6, .75, .83, .9, .98)
for(i in 1:9){
  text(vals[i], -.15, as.character(round(betas[i], 2)), cex = .8)
}

## @knitr gp-realizations

gr <- 100
xs <- seq(0, 1, length = gr)

# we'll use Matern correlation not exponential as exponential gives very wiggly sample paths
# note that with nu=0.5, the Matern is the same as the exponential correlation
maternCorr <- function(dist, rho, nu) {
    const <- -lgamma(nu)-(nu-1)*log(2)
    tmp <- exp(const+nu*log(2.0*sqrt(nu)*dist/rho)+log(besselK(2.0*sqrt(nu)*dist/rho,nu)))
    tmp[tmp>1] <- 1
    tmp[dist == 0] <- 1
    return(tmp)
}

mu <- 0
tau <- 1
rho <- c(.5, .7)
nu <- c(.5, 4)

library(fields)
dists <- rdist(xs)

# covariance matrices for two different sets of rho and nu
C <- L <- list()
for(i in 1:2) {
    C[[i]] <- tau^2 * maternCorr(dists, rho[i], nu[i])
    L[[i]] <- t(chol(C[[i]]))
}
    
set.seed(0)
par(mfrow  =  c(1, 2))

ncurves <- 5
z <- matrix(rnorm(gr*ncurves), ncol = ncurves)

plot(xs, mu + L[[1]] %*% z[ , 1], type = 'l', ylim = c(-2, 2),
     xlab = 'x', ylab = 'f(x)')
for(j in 2:ncurves)
    lines(xs, mu + L[[1]] %*% z[ , j], col = j)

plot(xs, mu + L[[2]] %*% z[ , 1], type = 'l', ylim = c(-2, 2),
     xlab = 'x', ylab = 'f(x)')
for(j in 2:ncurves)
    lines(xs, mu + L[[2]] %*% z[ , j], col = j)

## @knitr gp-fit

set.seed(0)
n <- 50
x <- seq(0, 1, len = n)
f <- sin(2*pi*x)

sigma <- 1
y <- rnorm(n, f, sigma)

plot(x, y)
lines(x, f)

library(fields)
library(nimble)

code <- nimbleCode({
    for(i in 1:n)
        y[i] ~ dnorm(f[i], sd = sigma)
    f[1:n] ~ dmnorm(mn[1:n], cov = C[1:n, 1:n])
    mn[1:n] <- mu*ones[1:n]
    sigma ~ dunif(0, 100)
    mu ~ dnorm(0, 0.0001)
    # exponential covariance probably induces functions that are too wiggly
    # but try it anyway
    C[1:n, 1:n] <- tau^2 * exp(-dists[1:n, 1:n] / rho)
    rho ~ dunif(0, 2)
    tau ~ dunif(0, 10)
})

m <- nimbleModel(code, constants = list(n = n,
                                        ones = rep(1, n),
                                        dists = rdist(x)),
                                        data = list(y = y),
                                        inits = list(mu = 0, rho = 0.5, tau = 1, sigma = 0.5))

cm <- compileNimble(m)

conf <- configureMCMC(m, monitors = c('rho', 'tau', 'sigma', 'mu', 'f'))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = m)

nIts <- 100000

cmcmc$run(nIts)

smp <- as.matrix(cmcmc$mvSamples)

postburn <- 40001:nIts
fCols <- grep("^f\\[.*\\]$", dimnames(smp)[[2]])
              
fMean <- apply(smp[postburn, fCols], 2, mean)

lines(x, fMean, col = 'red')
legend('topleft', lty = rep(1, 2), col = c('black','red'),
       legend = c('truth', 'fitted'))

                 


## @knitr finite-mixture-H

# let's look at number of important components using BDA-recommended (page 536)
# prior on mixture weights

library(nimble)  # nimble provides the rdirch function used here
m <- 1000
H <- 10
alpha0 <- 1
alpha <- rep(alpha0/H, H)
smp <- matrix(0, nrow = m, ncol = H)
for(i in 1:m)
    smp[i, ] <- rdirch(1, alpha)

par(mfrow = c(2,5))
for(h in 1:H)
    hist(smp[, h])

round(smp[1, ], 3)
round(smp[2, ], 3)


## @knitr setup-mixture

set.seed(0)
pi <- c(.4, .3, .2, .05, .03, .02)
realH <- length(pi)
mu <- rnorm(realH, 0, 5)
sigma <- runif(realH, 0.25, 1.5)

n <- 5000
c <- sample(1:realH, n, prob = pi, replace = TRUE)

mus <- mu[c]
sigmas <- sigma[c]

y <- rnorm(n, mus, sigmas)

doSmall <- FALSE
doTiny <- FALSE

if(doSmall) {
    n <- 500
    yFull <- y
    y <- y[1:n]
}
if(doTiny) {
    n <- 50
    yFull <- y
    y <- y[1:n]
}

# hist(y)


## @knitr fit-finite-mixture


deregisterDistributions('dmixN')

dmixN <- nimbleFunction(run = function(x = double(1),
                                       pi = double(1),
                                       mu = double(1), sigma = double(1),
                                       log = integer(0, default = 0)) {
    returnType(double(0))
    H <- length(pi)
    n <- length(x)
    logdens <- 0
    for(i in 1:n) {
        tmpdens <- 0
        for(h in 1:H) {
            tmpdens <- tmpdens + pi[h] * dnorm(x[i], mu[h], sigma[h], log = FALSE)
        }
        logdens <- logdens + log(tmpdens)
    }
        
    if(log) return(logdens) else return(exp(logdens))
})

rmixN <- nimbleFunction(run = function(n = integer(0), pi = double(1),
                       mu = double(1), sigma = double(1)) {
    stop("shouldn't be called")
    returnType(double(1))
    return(0)
})



code <- nimbleCode({
    y[1:n] ~ dmixN(pi[1:H], mu[1:H], sigma[1:H])
    # pi[1:H] ~ ddirch(alpha[1:H])
    for(h in 1:H) {
        alpha[h] <- alpha0 / H
        mu[h] ~ dnorm(mu0, sd = sigma0)
        sigma[h] ~ dlnorm(logmean, sdlog = logsigma0)
        piLatent[h] ~ dgamma(alpha[h], 1)
        pi[h] <- piLatent[h] / piLatentSum
    }
    piLatentSum <- sum(piLatent[1:H])
    mu0 ~ dnorm(0, sd = 3)
    sigma0 ~ dunif(0, 5)
    logmean ~ dnorm(0, sd = 2)
    logsigma0 ~ dunif(0, 2)
    alpha0 ~ dgamma(mean = 1, sd = 4)

})

set.seed(0)
Hmax <- 15
model <- nimbleModel(code, data = list(y = y),
                     constants = list(n = n, H = Hmax),
                     inits = list(piLatent = rep(1, Hmax),
                                  alpha0 = 1,
                                  mu = rnorm(Hmax, 0, 5),
                                  sigma = runif(Hmax, 0, 10),
                                  mu0 = 0,
                                  sigma0 = 3,
                                  logmean = 0,
                                  logsigma0 = 1
                                  ))

cmodel <- compileNimble(model)

mcmc <- buildMCMC(model, monitors = c('pi', 'mu', 'sigma', 'mu0',
                                      'sigma0', 'logmean', 'logsigma0', 'alpha0'))
cmcmc <- compileNimble(mcmc, project = model)

set.seed(0)
cmcmc$run(5000)



smp <- as.matrix(cmcmc$mvSamples)
save(smp, file = paste0('smpFinite', n, '.Rda'))

n <- 500
load(paste0('smpFinite', n, '.Rda'))
y <- yFull[1:n]

piCols <- grep("^pi.*\\]$", dimnames(smp)[[2]])
muCols <- grep("^mu.*\\]$", dimnames(smp)[[2]])
sigmaCols <- grep("^sigma.*\\]$", dimnames(smp)[[2]])


ts.plot(smp[, 'alpha0'])
par(mfrow = c(3,5))
for(i in piCols)
    ts.plot(smp[ , i])

par(mfrow = c(3,5))
for(i in muCols)
    ts.plot(smp[, i])

par(mfrow = c(3,5))
for(i in sigmaCols)
    ts.plot(smp[, i])


# look at density values on a grid of x values
gr <-  seq(-10, 10, length = 100)

dmixNcalc <- nimbleFunction(run = function(x = double(1),
                                       pi = double(2),
                                       mu = double(2), sigma = double(2),
                                       log = integer(0, default = 0)) {
    returnType(double(2))
    nIts <- dim(pi)[1]
    nGrid <- length(x)
    out <- matrix(nrow = nIts, ncol = nGrid)
    H <- dim(pi)[2]

    for(i in 1:nIts)
        for(j in 1:nGrid) {
            dens <- 0
            for(h in 1:H) 
                dens <- dens + pi[i, h] * dnorm(x[j], mu[i, h], sigma[i, h], log = FALSE)
            out[i, j] <- dens
        }
    return(out)
})

cdmixNcalc <- compileNimble(dmixNcalc)

smpDens <- cdmixNcalc(gr, smp[ , piCols], smp[ , muCols], smp[ , sigmaCols])

trueDens <- cdmixNcalc(gr, matrix(pi, nrow = 1), matrix(mu, nrow = 1),
                       matrix(sigma, nrow = 1))

par(mfrow = c(2,5))
subgrid <- seq(1, length(gr), by = 10)
for(j in subgrid)
    ts.plot(smpDens[ , j])

par(mfrow = c(1,1))
postBurn <- 501:nIts
hist(y, probability = TRUE, ncl = 40)
for(i in sample(postBurn, 10))
    lines(gr, smpDens[i, ], col = i)
lines(gr, colMeans(smpDens), lwd = 2)
lines(gr, trueDens, lwd = 2, lty = 2)

## @knitr fit-truncated-DPM

code <- nimbleCode({
    y[1:n] ~ dmixN(pi[1:H], mu[1:H], sigma[1:H])

    pi[1] <- V[1]
    prodterm[1] <- 1
    for(h in 2:(H-1)) {
        prodterm[h] <- prodterm[h-1] * (1-V[h-1]) 
        pi[h] <- V[h] * prodterm[h]
    }
    pi[H] <- 1 - sum(pi[1:(H-1)])
    for(h in 1:H) {
        V[h] ~ dbeta(1, alpha)
        mu[h] ~ dnorm(mu0, sd = sigma0)
        sigma[h] ~ dlnorm(logmean, sdlog = logsigma0)
    }
    mu0 ~ dnorm(0, sd = 3)
    sigma0 ~ dunif(0, 5)
    logmean ~ dnorm(0, sd = 2)
    logsigma0 ~ dunif(0, 2)
    
    alpha ~ dgamma(1, 0.25)
})

Hmax <- 15
alpha <- 2
set.seed(0)
V <- rbeta(Hmax, 1, alpha)

model <- nimbleModel(code, data = list(y = y),
                     constants = list(n = n, H = Hmax),
                     inits = list(V = V,
                                  alpha = 2,
                                  mu = rnorm(Hmax, 0, 5),
                                  sigma = runif(Hmax, 0, 10),
                                  mu0 = 0,
                                  sigma0 = 3,
                                  logmean = 0,
                                  logsigma0 = 1
                                  ))

cmodel <- compileNimble(model)

mcmc <- buildMCMC(model, monitors = c('pi', 'mu', 'sigma', 'alpha',
                                      'mu0', 'sigma0', 'logmean', 'logsigma0'))
cmcmc <- compileNimble(mcmc, project = model)

nIts <- 5000
set.seed(0)
cmcmc$run(nIts)

smp <- as.matrix(cmcmc$mvSamples)

save(smp, file = paste0('smpDPM', n, '.Rda'))

n <- 500
load(paste0('smpDPM', n, '.Rda'))
y <- yFull[1:n]

piCols <- grep("^pi.*\\]$", dimnames(smp)[[2]])
muCols <- grep("^mu.*\\]$", dimnames(smp)[[2]])
sigmaCols <- grep("^sigma.*\\]$", dimnames(smp)[[2]])

ts.plot(smp[ , 'alpha'])
par(mfrow = c(3,5))
for(i in piCols)
    ts.plot(smp[ , i])

par(mfrow = c(3,5))
for(i in muCols)
    ts.plot(smp[ , i])

par(mfrow = c(3,5))
for(i in sigmaCols)
    ts.plot(smp[ , i])

smpDens <- cdmixNcalc(gr, smp[ , piCols], smp[ , muCols], smp[ , sigmaCols])

trueDens <- cdmixNcalc(gr, matrix(pi, nrow = 1), matrix(mu, nrow = 1),
                       matrix(sigma, nrow = 1))

par(mfrow = c(2,5))
subgrid <- seq(1, length(gr), by = 10)
for(j in subgrid)
    ts.plot(smpDens[ , j])

par(mfrow = c(1,1))
hist(y, probability = TRUE, ncl = 40)
postBurn <- 501:nIts
for(i in sample(postBurn, 10))
    lines(gr, smpDens[i, ], col = i)
lines(gr, colMeans(smpDens), lwd = 2)      
lines(gr, trueDens, lwd = 2, lty = 2)
