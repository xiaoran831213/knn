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

source('sim/lnk.R')
test.lnk <- function()
{
    x <- rnorm(50000, 0, 1)
    dat <- within(list(),
    {
        NM <- DF(d='nn', x=x, y=x)
        st <- DF(d='st', x=x, y=dc(x, 'st', 0.01))
        ca <- DF(d='ca', x=x, y=dc(x, 'ca', 0.01))
        bn <- DF(d='bn', x=x, y=dc(x, 'bn', 0.00))
        ch <- DF(d='ch', x=x, y=dc(x, 'ch', 0.00))
        ps <- DF(d='ps', x=x, y=dc(x, 'ps', 0.00))
        ex <- DF(d='ex', x=x, y=dc(x, 'ex', 0.00))
    })
    dat <- do.call(rbind, dat)
    
    ## dst <- rbind(
    ##     DF(dst='nm', val=x),
    ##     DF(dst='st', val=x.st),
    ##     DF(dst='ca', val=x.ca))

    library(ggplot2)
    g <- ggplot()
    ## g <- g + geom_freqpoly(aes(color=dst), binwidth = 0.1)
    ## g <- g + xlim(-15, 15)
    ## g <- g + ylim(c(1, NA))
    g <- g + geom_line(aes(x=x, y=y, color=d), dat)
    g
}
