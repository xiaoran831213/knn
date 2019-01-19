lol2df <- function(x)
{
    dat <- lapply(list(...), function(l)
    {
        if(any(sapply(l, is.list)))
            l <- unlist(l)
        l
    })
    idx <- unique(unlist(sapply(dat, names)))
    dat <- lapply(dat, `[`, idx)
    dat
}

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

geom.mean <- function(x, na.rm=TRUE)
{
    (-1)^sum(x < 0, na.rm=na.rm) * exp(mean(log(abs(x)), na.rm=na.rm))
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
