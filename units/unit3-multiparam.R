## @knitr post-dep
library(MASS)

set.seed(0)
n <- 15
y <- rnorm(n, 0, 1)
ybar <- mean(y)
s2y <- var(y)

sig2 <- 1/rgamma(10000, (n-1)/2, (n-1)*s2y/2)
#  s2y*rchisq(10000, n-1)
mu <- rnorm(10000, ybar, sqrt(sig2/n))

par(mgp = c(1.8, .7, 0), mai = c(.6, .6, .4, .1), mfrow = c(1, 3))
pplot(sig2[1:1000], mu[1:1000], xlim = c(.4, 2.5), ylim = c(-1, 1), xlab = expression(sigma^2), ylab = expression(mu))
points(0, 1, pch = 'x', cex = 2)
kd <- kde2d(sig2, mu, n = 100)
image(kd, xlim = c(.4, 2.5), ylim = c(-1, 1), xlab = expression(sigma^2), ylab = expression(mu))
points(0, 1, pch = 'x', cex = 2)

fun <- function(param){
  mu <- param[2]
  sig2 <- param[1]
  alpha <- (n-1)/2; beta = (n-1)*s2y/2
  return(dnorm(mu, ybar, sqrt(sig2/n), log = TRUE)-(alpha+1)*log(sig2)-beta/sig2)
}

xs <- seq(.4, 2.5, len = 100)
ys <- seq(-1, 1, len = 100)
zs <- apply(expand.grid(xs, ys), 1, fun)
contour(xs, ys, matrix(zs, 100, 100), nlevels = 25, xlab = expression(sigma^2), ylab = expression(mu))

## @knitr gibbs

y <- c(78,66,65,63,60,60,58,56,52,50)
n <-  length(y)
ybar <- mean(y)
s2y <- var(y)

nIts <- 1000
mu <- sig2 <- rep(0, nIts)

# starting values
mu[1] <- ybar
sig2[1] <- s2y

set.seed(0)

for(t in 2:nIts){
   mu[t] <- rnorm(1, ybar, sqrt(sig2[t-1]/n))
   sig2[t] <- 1/rgamma(1, n/2, sum((y-mu[t])^2)/2)
}

muGrid <- seq(45, 75, length = 100)
sig2Grid <- seq(0, 300, length = 100)

dtgen <- function(x, mu, sigma, df = 1, log = FALSE) {
    result <- lgamma((df+1)/2) - lgamma(df/2) - 0.5*log(df*pi) - log(sigma)
    result <- result -((df+1)/2) * log(1 + ((x - mu)/sigma)^2 / df )
    if(log) return(result) else return(exp(result))
}

dinvgamma <- function(x, alpha, beta, log = FALSE) {
    result <- alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
    if(log) return(result) else return(exp(result))
}

#pdf('figures/gibbsEx.pdf')
par(mfrow = c(2,2))
hist(mu, main = expression(mu))
plot(density(mu))
lines(muGrid, dtgen(muGrid, ybar, sqrt(s2y/n), n-1), col = 'red')
hist(sig2, main = expression(sigma^2))
plot(density(sig2))
lines(sig2Grid, dinvgamma(sig2Grid, (n-1)/2, (n-1)*s2y/2), col = 'red')
#dev.off()
