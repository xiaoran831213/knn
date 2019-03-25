#' concatenate a list
CL <- function(ret=NULL, ...)
{
    dot <- list(...)
    ret <- ret %||% list()
    for(i in seq_along(dot))
    {
        ret <- c(ret, list(dot[[i]]))
    }
    idx <- seq(length(ret) - length(dot) + 1, l=length(dot))
    names(ret)[idx] <- names(dot)
    ret
}

#' extract elements from a list of lists
EL1 <- function(ll, keys, ret=c('list', 'data.frame'))
{
    ret <- match.arg(ret, c('list', 'data.frame'))
    ll <- lapply(ll, `[`, keys, drop=FALSE)
    ll <- lapply(ll, do.call, what=data.frame)
    if(ret == 'data.frame')
        ll <- do.call(rbind, ll)
    ll
}

EL2 <- function(ll, key, ret=c('list', 'data.frame'))
{
    ret <- match.arg(ret, c('list', 'data.frame'))
    ll <- lapply(ll, `[[`, key, drop=FALSE)
    if(ret == 'data.frame')
    {
        . <- lapply(ll, unlist)
        . <- lapply(., drop)
        if(all(sapply(sapply(., dim), is.null)))
            ll <- as.data.frame(do.call(rbind, .))
    }
    ll
}

mean.list <- function(x) Reduce(`+`, x) / length(x)

#' short hand for data.fram
DF <- data.frame

`%||%` <- function(x, y) if(is.null(x)) y else x

`%[%` <- function(ll, key) if(is.function(key)) lapply(ll, key) else lapply(ll, `[`, key)
`%$%` <- function(ll, key) if(is.function(key)) lapply(ll, key) else lapply(ll, `[[`, key)

#' vector row-binder
#'
#' row bind vectors of different elements
#' @param ... vectors or dafa.frames to be binded
rbd <- function(...)
{
    items <- list(...)
    items <- items[!sapply(items, is.null)]

    ## fix dimensions
    for(i in seq_along(items))
    {
        x <- items[[i]]
        r <- names(items)[[i]]
        if(length(dim(x)) != 2)
        {
            n <- names(x)
            dim(x) <- c(1L, length(x))
            colnames(x) <- n
        }
        if(nrow(x) == 1 && is.null(rownames(x)))
            rownames(x) <- r
        items[[i]] <- DF(x)
    }
    
    ## C=column names, and P=column counts
    C <- Reduce(union, sapply(items, colnames))
    P <- length(C)

    ## NA padding
    ret <- lapply(items, function(i)
    {
        n <- list(rownames(i), setdiff(C, colnames(i)))
        m <- matrix(NA, nrow(i), P - ncol(i), dimnames=n)
        cbind(i, m)[, C, drop=FALSE]
    })
    ret <- do.call(rbind, ret)
    ret
}

#' vector col-binder
#'
#' col bind vectors of different length
#' @param ... vectors to be binded
.cbd <- function(...)
{
    items <- list(...)
    nms <- Reduce(union, sapply(items, names))
    ret <- sapply(items, `[`, nms)
    rownames(ret) <- nms
    colnames(ret) <- names(items)
    ret
}

#' Dimension Change.
.dim <- function(x, ...) {dim(x) <- c(...); x}

#' Multiplied by Matrix.
.mbm <- .Primitive('%*%')

#' Inverse Symmetric Matrix.
.inv <- function(x, ...) chol2inv(chol(x))

#' format print to STDOUT
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
