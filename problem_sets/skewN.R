dskewN <- nimbleFunction(
    # see Arellano-Valle and Azzalini; The centred parametrization for the multivariate skew-normal distribution; Journal of Multivariate Analysis 99 (2008) 1362–1382
    run = function(x = double(0), mu = double(0), sigma = double(0),
                   gamma1 = double(0),
        log = integer(0, default = 0)) {
        returnType(double(0))
        
        pi <- 3.141593
        pos <- gamma1 > 0
        gamma1 <- abs(gamma1)
        c <- (2*gamma1 / (4-pi)) ^ (1/3)
        muz <- c / sqrt(1+c^2)
        delta <- muz / sqrt(2/pi)
        alpha <- delta / sqrt(1-delta^2)
        
        if(!pos) {
            alpha <- -alpha
            delta <- -delta
            muz <- -muz
        }
        
        omega <- sigma/sqrt(1 - muz^2)
        xi <- mu - omega*muz
        
        z <- (x-xi)/omega
        logProb <- log(2) - log(omega) + dnorm(z, 0, 1, log = TRUE) +
            pnorm(alpha * z, 0, 1, log.p = TRUE)
        if(log) return(logProb)
        else return(exp(logProb))
    })

# we need this if we want to simulate from the prior
# code based on rsn() from sn package
rskewN  <- nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0),
                   gamma1 = double(0)) {
        returnType(double(0))

        pi = 3.141593
        sign <- gamma1 > 0
        gamma1 <- abs(gamma1)
        
        c = (2*gamma1 / (4-pi)) ^ (1/3)
        muz = c / sqrt(1+c^2)
        
        delta = muz / sqrt(2/pi)
        
        if(!sign)  {
            delta <- -delta
            muz <- -muz
        }
        
        omega <- sigma/sqrt(1 - muz^2)
        xi <- mu - omega*muz
        
        truncN <- qnorm(runif(1, 0.5, 1))
        z <- delta * truncN + sqrt(1 - delta^2) * rnorm(1)
         
        return(xi + omega * z)
    })

registerDistributions(list(
        dskewN = list(
               BUGSdist = "dskewN(mu, sigma, gamma1)")))
