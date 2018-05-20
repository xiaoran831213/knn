cl <- function(ret=NULL, ...)
{
    dots <- list(...)
    ret <- ret %||% list()
    for(i in seq_along(dots))
    {
        ret <- c(ret, list(dots[[i]]))
    }
    ret
}
df <- data.frame
DF <- data.frame

`%||%` <- function(x, y) if(is.null(x)) y else x

#' Dimension Change.
.dim <- function(x, ...) {dim(x) <- c(...); x}

#' Multiplied by Matrix.
.mbm <- .Primitive('%*%')

#' Inverse Symmetric Matrix.
.inv <- function(x, ...) chol2inv(chol(x))

PF <- function(...) cat(sprintf(...))

str.cfg <- function(cfg)
{
    cfg <- unlist(cfg)
    paste0(names(cfg), '=', cfg, collapse=', ')
}

nlk <- function(x, v=NULL, u=NULL, a=NULL, ...)
{
    v <- v %||% diag(length(x))
    u <- u %||% chol(v)
    a <- a %||% chol2inv(u)
    .5 * sum(crossprod(x, a) * x) + sum(log(diag(u))) + .5 * nrow(v) * log(2*pi)
}

lmm <- function(y, v, e)
{
    N <- nrow(v)
    f <- v - diag(e, N)
    u <- chol(v)
    a <- chol2inv(u)

    ## prediction 1: conditional Gaussian
    h <- f %*% (a %*% y)
    mse <- mean((y - h)^2)
    cyh <- cor(y, h)

    ## prediction 2: leave one out CV
    ## h <- y - a %*% y / diag(a)
    ## loo <- mean((y - h)^2)

    ## negative log likelihood
    nlk <- nlk(y, v, u, a)

    df(key=c('mse', 'nlk', 'cyh'), val=c(mse, nlk, cyh))
}

nul <- function(y)
{
    N <- NROW(y)

    ## the distribution of y
    m <- mean(y)
    v <- diag(var(y), N)

    ## mean square error
    mse <- mean((y - m)^2)

    ## negative log likelihood
    nlk <- nlk(y, v)

    ## leave one out
    ## h <- (sum(y) - y) / (N - 1)
    ## loo <- mean((y - h)^2)

    ## report
    df(key=c('mse', 'nlk'), val=c(mse, nlk))
}


#' @title multi variant norm from upper Cholesky
#' Most mathimatical notation uses lower Cholesky Decomposition (CD), but
#' upper CD is faster for most computer programs.
#' MVN(p) := xu, where x ~ N(0, I(p)), {p} is the dimension of Sigma, and
#' {u} is the upper CD of Sigma.
#' @param n the number of multivariant norm samples to draw;
#' @param v the covariance matrix, or its cholesky decomposition.
#' @param u TRUE to treat \code{v} as a cholescky decomposition.
#' @return MVN samples organized in an N-P matrix.
mvn <- function(n, v, u=NULL, ...)
{
    u <- if(u %||% FALSE) v else chol(v)
    p <- nrow(u)
    x <- matrix(rnorm(n * p), n, p) %*% u
    x
}

xvy <- function(x, v, y=x)
{
    sum(crossprod(x, v) * y)
}

softplus <- function(x) log(1+exp(x))
logistic <- function(x) 1/(1+exp(-x))

softmax <- function(x)
{
    e <- exp(x - max(x))
    e / sum(e)
}

softxam <- function(x, p=NULL)
{
    if(is.null(p))
        p <- softmax(x)
    p * (diag(length(x)) - p)
}

looCV <- function(x, cv, ac=NULL, va=FALSE, ...)
{
    if(is.null(ac))
    {
        dc <- chol(cv)
        ac <- chol2inv(dc)
    }
    h <- x - ac %*% x / diag(ac)
    if(va)
        ret <- list(h=h, v=1/diag(ac))
    else
        ret <- h
    ret
}


## generate random rotation
rrot <- function(p=2, m=0, s=1)
{
    a <- c(1, rep(0, p - 1))
    b <- a + rnorm(p, m, s)
    b <- b / sqrt(sum(b^2))
    ab <- sum(a * b)
    ca <- a - b * ab
    ca <- ca/sqrt(sum(ca^2))
    A <- b %*% t(ca)
    A <- A - t(A)
    theta <- acos(ab)
    diag(p) + sin(theta) * A + (cos(theta) - 1) * (b %*% t(b) + ca %*% t(ca))
}
