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
    cdc <- as.logical(get0('cdc', ifnotfound=TRUE))
    print(list(cdc=cdc))
    if(cdc)
        z <- z - y
    z
}

NL <- function(x) x
SC <- function(x) drop(scale(x, scale=FALSE))
ST <- function(x) dc(x, 'st', 0.05)
CA <- function(x) dc(x, 'ca', 0.05)
BN <- function(x) dc(x, 'bn')
PS <- function(x) dc(x, 'ps')
XP <- function(x) dc(x, 'ex')
X2 <- function(x) dc(x, 'ch')

P1 <- function(x) (1 + SC(x))^1
P2 <- function(x) (1 + SC(x))^2
P3 <- function(x) (1 + SC(x))^3
O1 <- function(x) SC(x)^1
O2 <- function(x) SC(x)^2
O3 <- function(x) SC(x)^3

SN <- function(x) sin(2 * pi * x)
