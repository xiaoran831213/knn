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

test.out3 <- function(N=500, P=1000, t=20)
{
    library(microbenchmark)
    a <- matrix(rnorm(N * P), N, P)

    f1 <- function()
    {
        (apply(a, 1L, function(x) sum(tcrossprod(x, x))) - rowSums(a^2))/2
    }
    
    f2 <- function()
    {
        r <- numeric(N)
        for(i in seq.int(P-1))
        {
            for(j in seq.int(i+1, P))
            {
                r <- r + a[, i] * a[, j]
            }
        }
        r
    }

    print(microbenchmark(r1 <- f1(), r2 <- f2(), times=t))
    
    list(r1, r2)
}

test.kpl1 <- function(N=500, P=2000)
{
    g <- readRDS('data/p35_c05.rds')$gmx
    x <- g[sample.int(nrow(g), N), sample.int(ncol(g), P)]
    x
}
