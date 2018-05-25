## test kenrl minque
library(MASS)
library(microbenchmark)
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")

test <- function(N=100, P=200, r=100)
{
    x <- matrix(rnorm(N * P), N, P)
    v <- list(e=diag(nrow(x)), p=tcrossprod(x))
    k <- length(v)
    p <- rbind(diag(k), rep(1, k))

    svc <- .1 * v$e + 1 * v$p
    y <- mvrnorm(1, rep(0, N), svc)

    bmk <- microbenchmark(
        r1 <- .Call('_knn_knl_mnq_cpp', PACKAGE = 'knn', as.matrix(y), v, p),
        r2 <- knl.mnq(y, v),
        times = r)
    print(bmk)
    invisible(list(r1=r1, r2=r2))
}
