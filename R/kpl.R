## main function

#' kernel executor
#'
#' @param x N-P matrix or N-1 vector, the data on which to apply the
#' kernel functions
#' @param k a kernel function or list of kernel functions
#' @param ... additional named parameters for the kernel function(s)
#' @return a list of N-N kernel matrices
krn <- function(x, M=~ply, J=FALSE)
{
    N <- NROW(x)
    enlist <- function(x) if(is.list(x)) x else list(x)

    ## function to recusively retrive upper triangles
    mtx <- matrix(0, N, N)
    utr <- upper.tri(mtx, TRUE)
    upt <- function(x) if(is.list(x)) lapply(x, upt) else x[utr]

    ## get kernel functions
    tm <- terms(M)
    vs <- eval(attr(tm, 'variables'))
    vs <- do.call(c, vs)
    names(vs) <- rownames(attr(tm, 'factors'))
    
    ## evaluate basic kernels
    ks <- lapply(vs, function(v) v(x))
    
    ## get upper triangles
    ks <- upt(ks)
    
    ## expand formula
    pt <- gsub('([().])', '[\\1]', names(ks))
    nm <- lapply(ks, names)
    rp <- sapply(nm, function(n) sprintf('(%s)', paste0(n, collapse='+')))
    
    M <- as.character(M)
    for(i in seq_along(ks))
        M <- sub(pt[i], rp[i], M)
    M <- as.formula(M)
    
    ## expand data
    ks <- do.call(c, ks)
    names(ks) <- do.call(c, nm)
    ks <- model.matrix(M, ks)

    ## include J kernel?
    if(J) colnames(ks)[1] <- "J" else ks <- ks[, -1, drop=FALSE]
    
    ## reconstruct kernel matrices
    ks <- as.data.frame(ks)
    lapply(ks, function(k) {mtx[utr] <- k; mtx <- t(mtx); mtx[utr] <- k; mtx})
}

knl2str <- function(k)
{
    if(!is.list(k))
        k <- lapply(k, as.list)
    ## covert all kernels to list of lists
    if(any(sapply(k, is.list)))         # at least one kernel alreay in list format
        k <- lapply(k, as.list)
    else                                # a single kernel in list format
        k <- list(k)
    
    k <- lapply(k, function(.)
    {
        . <- unlist(.)
        p <- .[-1]
        p <- paste0(paste(names(p), p, sep='='), collapse=', ')
        p <- paste0('(', p, ')')
        sprintf('%s%s', .[1], p)
    })
    k <- paste(k, collapse=' + ')
    k
}

#' identity kernel function in matrix form
#'
#' @param x N-P matrix or N-1 vector of data
#' @param y not used
#' @param ... not used
#' 
#' @return an N-N identity matrix
idn <- function(x, y=NULL, ...)
{
    diag(NROW(x))
}
    
#' gaussian kernel function in matrix form
#'
#' @param x N-P matrix or N-1 vector of data
#' @param y not used
#' @param s variation length scale
#' @param g scaling factor
#' @param ... not used
#' 
#' @return an N-N identity matrix
gau <- function(x, y=NULL, s=1, g=1/NCOL(x), ...)
{
    exp(-euc2(x, y) * (.5 * g / s^2))
}

lap <- function(x, y=NULL, s=1, g=1/NCOL(x), ...)
{
    exp(-as.matrix(dist(x, 'man')) * g / s)
}


#' @title polynomial kernel between X and Y
#' @param x: matrix N rows of samples, P columns of features
#' @param y: matrix M rows of samples, P columns of features
#' @param g: numberic, default is 1/P
#' @param coef0: numberic, default is 0
#' @param degree: integer, default is 1
#' @details K(X, Y) = (g <X, Y> + coef0)^degree;
#' by default, coef=0, degree=1, g=1/P;
#' when Y is is NULL, K(X, X) is computed.
ply <- function(x, y=NULL, g=1/NCOL(x), coef0=0, degree=1L, ...)
{
    (tcrossprod(x, y) * g + coef0) ^ degree
}

#' @title polynomial multi-way interaction
#' @param x: matrix N rows of smaples, P columns of features
#' @param y: matrix M rows of smaples, P columns of features
pqw <- function(x, y=NULL, g=1/NCOL(x), q=2, ...)
{
    if(is.null(y))
        oq <- tcrossprod(x^q)
    else
        oq <- tcrossprod(x^q, x^q)
    (g^q) * (tcrossprod(x, y)^q - oq)
}

p2w <- function(x, y=NULL, g=1/NCOL(x), ...)
{
    if(is.null(y))
        o2 <- tcrossprod(x^2)
    else
        o2 <- tcrossprod(x^2, y^2)
    (g^2) * (tcrossprod(x, y)^2 - o2)
}

p2o <- function(x, y=NULL, g=1/NCOL(x), ...)
{
    if(is.null(y))
        (g^2) * tcrossprod(x^2)
    else
        (g^2) * tcrossprod(x^2, y^2)
}

#' Identity by State Kernel
#'
#' For genomic dosage data, the similarity between two individual
#' i and j, contributed to the k th. variants is:
#'
#' $ s_{ij} = 2 - |g_{ik} - g_{jk}| $
#'
#' @param x N matrix of N row individuals and P column variants
#' @param l level of genomic variation (def = 2, i.e., allele dosage.)
#' @return NxN IBS kernel matrix
ibs <- function(x, l=2, ...)
{
    if(is.null(l))
        x <- apply(x, function(.) . / max(.))
    else
        x <- x / l
    1 - as.matrix(dist(x, 'man')) / ncol(x)
}

#' Genetic Relatedness Matrix
#'
#' The kernel is enssentially a product kernel weighted
#' by estimated standard deviation on each genomic variant,
#' as Yang et.al stated in their 2014 paper for GCTA.
#' 
#' @param x matrix of N row individuals and P colunn variant
#' @param o integer order of polynomial order of the kernel
#' @param q vector of user provided weights, if null, the
#' weights are estimated standard deviation.
#'
#' @return NxN kernel matrix
grm <- function(x, o=1, q=NULL)
{
    if(is.null(q))                      # AF
        q <- colMeans(x, TRUE) / 2
    s <- sqrt((2 * q * (1 - q)))^o      # weights

    x <- x[, s > 0]                     # remove degeneracy
    q <- q[  s > 0]                     #
    s <- s[  s > 0]                     # 
    
    a <- is.na(x)
    M <- tcrossprod(1 - a)              # pairwise non-NA
    x <- as.matrix(scale(x, q * 2, s))

    ## set NA to zero
    x[a] <- 0.0
    rm(a)

    k <- tcrossprod(x) / M
    ## k <- k / mean(diag(k))
    k
}

cmb <- function(k, w, drop=FALSE)
{
    w <- as.matrix(w)
    stopifnot(is.list(k), length(k) == NROW(w))

    M <- ncol(w)                        # num of output
    L <- nrow(w)                        # num of kernel
    m <- 1
    r <- list()
    for(m in seq.int(1, M))
    {
        . <- k[[1]] * w[1, m]
        l <- 2
        while(l <= L)
        {
            . <- . + k[[l]] * w[l, m]
            l <- l + 1
        }
        r[[m]] <- .
        m <- m + 1
    }

    ## drop the list if only one output is produced?
    if(drop && M==1)
        r <- r[[1]]
    r
}

## kinship matrix
kin <- function(x, y=NULL)
{
    N <- NROW(x)                       # sample size
    ## -1 for AA, 0 for Aa, 1 for aa
    x <- x - 1
    if(!is.null(y))
        y <- y - 1
    k <- tcrossprod(x, y)
    k / sum(diag(k)) * N
}

ntk <- function(x, y=NULL, s0=1, s1=1)
{
    x <- cbind(1, x)
    if(is.null(y))
        y <- x
    else
        y <- cbind(1, y)
    P <- ncol(x)
    N <- nrow(x)
    S <- diag(c(s0, rep(s1, P-1)))

    xx <- x %*% S %*% t(x)
    yy <- y %*% S %*% t(y)
    xy <- x %*% S %*% t(y)

    z <- 2 * xy / sqrt((1 + 2 * xx) * (1 + 2 * yy))
    2 / pi * asin(z)
}

ack <- function(x, y=NULL)
{
    x2 <- rowSums(x^2)
    y2 <- if(is.null(y)) x2 else rowSums(y^2)
    xy <- tcrossprod(x, y)
    n2 <- outer(x2, y2)
    acos(xy/sqrt(outer(x2, y2)))
}

psd <- function(K)
{
    with(eigen((K + t(K)) / 2), vectors %*% (pmax(values, 0) * t(vectors)))
}

std <- function(K)
{
    K / sum(diag(K)) * NROW(K)
}
