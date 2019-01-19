library(plinkBED)
kpa <- function(x)
{
    x <- na.omit(x)
    i <- rownames(x)
    x <- x[, i]
    kappa(x, method="direct")
}


swt <- function(x, o=1, q=NULL)
{
    if(is.null(q))                      # AF
        q <- colMeans(x, TRUE) / 2
    s <- sqrt((2 * q * (1 - q)))^o      # weights


    x <- x[, s > 0]                     # remove degeneracy
    q <- q[  s > 0]                     #
    s <- s[  s > 0]                     # 
    
    a <- is.na(x)
    M <- tcrossprod(1 - a)              # pairwise non-NA
    x <- as.matrix(scale(x, q * 2, s))

    ## set NA to zero
    x[a] <- 0.0
    rm(a)

    k <- tcrossprod(x) / M
    ## k <- k / mean(diag(k))
    k
}

S1T <- function(x) swt(x, 1)
S2T <- function(x) swt(x, 2)
S3T <- function(x) swt(x, 3)
S4T <- function(x) swt(x, 4)
