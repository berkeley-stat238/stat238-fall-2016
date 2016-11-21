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

## @knitr gp-example

gr <- 100
xs <- seq(0, 1, length = gr)

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


## @knitr finite-mixture-H

library(nimble)
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



## @knitr fit-dpm

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

# hist(y)


## @knitr fit-finite-mixture

deregisterDistributions('dmixN')

dmixN <- nimbleFunction(run = function(x = double(1),
                                       pi = double(1),
                                       mu = double(1), sigma = double(1), log = integer(0, default = 0)) {
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
    pi[1:H] ~ ddirch(alpha[1:H])
    for(h in 1:H) {
        alpha[h] <- alpha0 / H
        mu[h] ~ dnorm(mu0, sd = sigma0)
        sigma[h] ~ dlnorm(logmean, sdlog = logsigma0)
    }
    mu0 ~ dnorm(0, sd = 3)
    sigma0 ~ dunif(0, 5)
    logmean ~ dnorm(0, sd = 2)
    logsigma0 ~ dunif(0, 2)
    alpha0 ~ dgamma(mean = 1, sd = 4)
})

Hmax <- 25
model <- nimbleModel(code, data = list(y = y),
                     constants = list(n = n, H = Hmax),
                     inits = list(pi = rep(1/Hmax, Hmax),
                                  alpha0 = 1,
                                  mu = rnorm(Hmax, 0, 5),
                                  sigma = runif(Hmax, 0, 10),
                                  mu0 = 0,
                                  sigma0 = 3,
                                  logmean = 0,
                                  logsigma0 = 1
                                  ))

cmodel <- compileNimble(model)

mcmc <- buildMCMC(model)
cmcmc <- compileNimble(mcmc, project = model)

cmcmc$run(1000)

smp <- as.matrix(cmcmc$mvSamples)

                               
