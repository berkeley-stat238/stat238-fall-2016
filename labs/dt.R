dt <- function(x, mu, sigma, df = 1, log = FALSE) {
    result <- lgamma((df+1)/2) - lgamma(df/2) - 0.5*log(df*pi) - log(sigma)
    result <- result -((df+1)/2) * log(1 + ((x - mu)/sigma)^2 / df )
    if(log) return(result) else return(exp(result))
}
