library(plinkBED)
kpa <- function(x)
{
    x <- na.omit(x)
    i <- rownames(x)
    x <- x[, i]
    kappa(x, method="direct")
}

S1T <- function(x) swt(x, 1)
S2T <- function(x) swt(x, 2)
S3T <- function(x) swt(x, 3)
S4T <- function(x) swt(x, 4)
