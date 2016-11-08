## @knitr multiple-testing

m <- 10000
y <- rnorm(m, 0, 1)

pvalue <- 2*pnorm(-abs(y))
sum(pvalue < 0.05)

# Benjamini-Hochberg

qstar = 0.05
sorted_ps <- sort(pvalue)
which(sorted_ps <= (1:m)*qstar / m)

# Bayesian formulation

# suppose abs(theta) > 1.5 is considered 'interesting'

scientificSign <- 1.5

# basic Gibbs

S <- 5000

sigma2 <- 1
tau2 <- rep(0, S)
theta <- matrix(0, S, m)
mu <- rep(0, S)

tau2[1] <- 1
theta[1, ] <- rnorm(m, 0, sqrt(tau2[1]))

for(s in 2:S) {
    tau2[s] <- 1/rgamma(1, -0.5 + m/2, sum((theta[s-1,]-mu[s-1])^2) /2 )
    mu[s] <- rnorm(1, mean(theta[s-1,]), sqrt(tau2[s]/m))
    vInv <- (1/tau2[s] + 1/sigma2)
    theta[s, ] <- rnorm(m, (mu[s]/tau2[s] + y/sigma2)/vInv, sqrt(1/vInv))
}

postp <- apply(theta, 2, function(x) mean(abs(x) > scientificSign))
thetaMns <- colMeans(theta)
thetaSds <- apply(theta, 2, sd)


thetaNew <- matrix(0, S, m)
mu <- 0

tau2 <- 10000
thetaNew[1, ] <- rnorm(m, 0, 1)

for(s in 2:S) {
    vInv <- (1/tau2 + 1/sigma2)
    thetaNew[s, ] <- rnorm(m, (mu/tau2 + y/sigma2)/vInv, sqrt(1/vInv))
}

# "scientifically significant"
postp <- apply(thetaNew, 2, function(x) mean(abs(x) > scientificSign))
sum(postp > 0.9)
sum(postp > 0.95)

# Bayesian quasi-p-value
# find thetas where there is substantial posterior probability they are positive or negative
ppos <- apply(thetaNew, 2, function(x) mean(x > 0))
pneg <- apply(thetaNew, 2, function(x) mean(x < 0))
pvalue_bayes <- apply(cbind(ppos, pneg), 1, max)

sum(pvalue_bayes > .95)


