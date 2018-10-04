## link functions
#' Distribution Convertor
#'
#' covert an assumed empirical nromal distribution to another distribution.
#'
#' x: the data that descript the empirical normal
#' d: the target distribution
#' .: additional paramters required by the target quantile (e.g., degree
#' of freedom for t and chisq, min and min for unif, etc.)
dc <- function(x, d=NULL, curb=0.01, ...)
{
    v <- var(x)                         # variance
    s <- sd(x)                          # sd
    m <- mean(x)                        # mean
    n <- length(x)                      # numbers
    p <- pnorm(x, m, s)                 # quantile

    p <- curb / 2 + p * (1 - curb)

    ## mean of target distribution.
    y <- switch(d, bn={2 * v}, ps={v}, ex={1 / s}, ch={v}, 0)

    z <- switch(
        d,
        ## st={     qt(p, (2 * v) / max(v - 1, 0))},
        st={     qt(p, 1 / s)},
        bn={ qbinom(p, ceiling(4 * v), .5)},
        ps={  qpois(p, v)},
        ex={   qexp(p, 1 / s)},
        ## ca={qcauchy(p) -> .; s * . / sd(.)},
        ca={qcauchy(p, 0, s)},
        ch={ qchisq(p, v / 2)})

    ## centered?
    cdc <- as.logical(get0('cdc', ifnotfound=FALSE))
    print(list(cdc=cdc))
    if(cdc)
        z <- z - y
    z
}

ST <- function(x) dc(x, 'st', 0.05)
BN <- function(x) dc(x, 'bn')
PS <- function(x) dc(x, 'ps')
EX <- function(x) dc(x, 'ex')
CA <- function(x) dc(x, 'ca', 0.05)
CH <- function(x) dc(x, 'ch')
NL <- function(x) x
I2 <- function(x) scale(drop((x + 1)^2))
I3 <- function(x) scale(drop((x + 1)^3))
O1 <- function(x) x
O2 <- function(x) drop(scale(x^2))
O3 <- function(x) drop(scale(x^3))
SN <- function(x) sin(2 * pi * x)
