#' Normal Distribution Convertor
#'
#' covert an assumed empirical nromal distribution to another distribution.
#'
#' x: the data that descript the empirical normal
#' q: the quantile function of the target distribution
#' .: additional paramters required by the target quantile (e.g., degree
#' of freedom for t and chisq, min and min for unif, etc.)
ndc <- function(x, q=qnorm, ...) q(pnorm(x, mean(x), sd(x)), ...)

## fast (squared) Euclidean distance
euc2 <- function(x, y=NULL)
{
    x2 <- rowSums(x^2)
    y2 <- if(is.null(y)) x2 else rowSums(y^2)
    xy <- tcrossprod(x, y)

    outer(x2, y2, `+`) - 2 * xy
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

