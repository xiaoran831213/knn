#' get a consecutive section from the genomic matrix
#' for training and testing dataset.
#' 
#' @param the basic gnomic matrix
#' @param N the size of development data
#' @param H the size of evaluation data
#' @param P the number of variants
get.gmx <- function(gmx, N, H=N, P=ncol(gmx) * .5)
{
    H <- min(N, nrow(gmx) - N)
    P <- min(P, ncol(gmx))
    
    ## individual indices
    idx <- sample.int(nrow(gmx), N + H)
    ## variant indices
    jdx <- seq(sample.int(ncol(gmx) - P, 1), l=P)

    ## developing data
    gmx.dvp <- as.matrix(gmx[idx[+(1:N)], jdx])

    ## evaluating data
    gmx.evl <- as.matrix(gmx[idx[-(1:N)], jdx])
    
    list(dvp=gmx.dvp, evl=gmx.evl)
}

## Hutchinson trace for symmetric matrix product
tr.hut <- function(..., N=100)
{
    lmx <- list(...)
    s <- c(-1, +1)
    n <- nrow(lmx[[1]])
    a <- 0
    for(i in 1:N)
    {
        z <- sample(s, n, TRUE)
        a <- a + sum(Reduce(tcrossprod, lmx, z) * z)
    }
    a / N
}

tr <- function(...)
{
    lmx <- list(...)
    prd <- Reduce(`%*%`, lmx)
    sum(diag(prd))
}

## condjugated gradient decent solver of Kx=b
sv.cgd <- function(K, b, max.itr=1000, tol=1e-6)
{
    s0 <- rnorm(length(b), sd=.1)
    dim(s0) <- dim(b)
    e0 <- b - K %*% s0
    d0 <- e0
    m0 <- sum(e0^2)
    for(i in seq.int(max.itr))
    {
        kd <- K %*% d0
        a0 <- m0 / sum(d0 * kd)
        s1 <- s0 + a0 * d0
        e1 <- e0 - a0 * kd
        m1 <- sum(e1^2)
        if(mean(m1)/N < tol)
            break
        s0 <- s1
        e0 <- e1
        d0 <- e1 + m1/m0 * d0
        m0 <- m1
    }
    as.numeric(s1)
}


test.solve <- function(N=1000)
{
    K <- matrix(rnorm(N^2), N, N)
    K <- tcrossprod(K)
    b <- rnorm(N)

    mb <- microbenchmark(
        r1 <- solve(K, b),
        r2 <- cg.solve(K, b, 1e5), times=2)
    print(mb)

    sum((r1 - r2)^2)
}
