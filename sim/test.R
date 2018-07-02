test.euc2 <- function(N=500, P=1000, t=20)
{
    library(microbenchmark)

    a <- matrix(rnorm(N * P), N, P)
    r <- microbenchmark(
        d1 <- as.matrix(dist(a))^2,
        d2 <- sqrt(abs(euc2(a)))^2,
        times=t)
    print(r)
    
    sum(abs(round(d1 - d2, 6)))
}


test.solve <- function(N=500, P=1000, t=20)
{
    library(microbenchmark)
    a <- tcrossprod(matrix(rnorm(N * P), N, P))
    b <- rnorm(N)

    inv <- function() {solve(a) %*% b}
    slv <- function() {solve(a, b)}
    qrs <- function() {qr.solve(a)}
    chl <- function() {chol2inv(chol(a)) %*% b}
    kpa <- function() {kappa(a)}

    r <- microbenchmark(inv, qrs, slv, chl, kpa, times=t)
    r
}
