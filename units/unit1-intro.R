## @knitr pri_to_post

n <- 15 
y<- 6
theta <- seq(0, 1, len = 100) 

# first example: uniform prior
alpha <- beta <- 1
pri1<- dbeta(theta, alpha, beta)
post1 <- dbeta(theta, alpha+y, beta+n-y) 

# second example: informative prior
alpha <- 2; beta <- 7
post2 <- dbeta(theta, alpha+y, beta+n-y)
pri2 <- dbeta(theta, alpha, beta)

# first example
par(mfrow = c(1, 2))
plot(theta, 20*dbinom(y, n, theta), xlab=expression(theta),  ylab='density', type='l', ylim=c(0, 5),  main = 'uniform prior case') 
lines(theta, post1, lty=1, lwd=2)
lines(theta, pri1, lty=1, col='darkgrey')
legend("topleft",  legend = c('prior',  'lik',  'posterior'), lty = rep(1, 3), col = c('darkgrey', 'black','black'), lwd = c(1, 1, 2), bty="n") 

# second example
plot(theta, 20*dbinom(y, n, theta), xlab=expression(theta), ylab='',type='l', ylim=c(0, 5),  main = 'informative prior case')
lines(theta, post2, lwd=2) 
lines(theta, pri2, col='darkgrey') 
legend("topleft", legend = c('prior', 'lik', 'posterior'), lty = rep(1,3), col = c('darkgrey', 'black','black'), lwd = c(1,1,2), bty="n") 
