## @knitr compare_priors

alpha1 <- beta1 <- 1

alpha2 <- beta2 <- 0

alpha3 <- beta3 <- 1/2

theta <- seq(0, 1, len = 100) 

# priors

pri1 <- dbeta(theta, alpha1, beta1) # unif orig scale
pri2 <- 0.1 * 1/(theta*(1-theta))   # unif logit scale
pri3 <- dbeta(theta, alpha3, beta3) # jeffreys

par(mfrow = c(1,3))
plot(theta, pri1, xlab=expression(theta),  ylab='density',
     type = 'l', ylim = c(0,5), main = 'priors')
lines(theta, pri2, lty = 2)
lines(theta, pri3, lty = 3)

legend("topleft", legend = c("uniform-prob", "uniform-logit", "jeffreys"), lty=1:3, bty='n')

# limited data

n <- 6
y<- 2

post1 <- dbeta(theta, alpha1+y, beta1 + n-y)
post2 <- dbeta(theta, alpha2+y, beta2 + n-y)
post3 <- dbeta(theta, alpha3+y, beta3 + n-y)

plot(theta, post1, xlab=expression(theta),  ylab='density',
     type = 'l', ylim = c(0,5), main = 'posteriors, small n')
lines(theta, post2, lty = 2)
lines(theta, post3, lty = 3)

legend("topleft", legend = c("uniform-prob", "uniform-logit", "jeffreys"), lty=1:3, bty='n')

# rich data

n <- 60
y<- 20

post1 <- dbeta(theta, alpha1+y, beta1 + n-y)
post2 <- dbeta(theta, alpha2+y, beta2 + n-y)
post3 <- dbeta(theta, alpha3+y, beta3 + n-y)

plot(theta, post1, xlab=expression(theta),  ylab='density',
     type = 'l', ylim = c(0,8), main = 'posteriors, larger n')
lines(theta, post2, lty = 2)
lines(theta, post3, lty = 3)

legend("topleft", legend = c("uniform", "uniform-logit", "jeffreys"), lty=1:3, bty='n')

