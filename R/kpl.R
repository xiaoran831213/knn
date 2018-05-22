## main function

#' kernel executor
#'
#' @param x N-P matrix or N-1 vector, the data on which to apply the
#' kernel functions
#' @param k a kernel function or list of kernel functions
#' @param ... additional named parameters for the kernel function(s)
#' @return a list of N-N kernel matrices
krn <- function(x, k=c(idn, ply), ...) lapply(k, do.call, list(x))

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
#' @param sigma variation length scale
#' @param gamma scaling factor
#' @param ... not used
#' 
#' @return an N-N identity matrix
gau <- function(x, y=NULL, sigma=.1, gamma=1/NCOL(x), ...)
{
    exp(-as.matrix(dist(x, 'euc') * gamma) / (2 * sigma^2))
}

lap <- function(x, y=NULL, sigma=.1, gamma=1/NCOL(x), ...)
{
    exp(-as.matrix(dist(x, 'man') * gamma) / (1 * sigma^1))
}


#' @title polynomial kernel between X and Y
#' @param x: matrix N rows of samples, P columns of features
#' @param y: matrix M rows of samples, P columns of features
#' @param gamma: numberic, default is 1/P
#' @param coef0: numberic, default is 0
#' @param degree: integer, default is 1
#' @details K(X, Y) = (gamma <X, Y> + coef0)^degree;
#' by default, coef=0, degree=1, gamma=1/P;
#' when Y is is NULL, K(X, X) is computed.
ply <- function(x, y=NULL, gamma=1/NCOL(x), coef0=0, degree=1L, ...)
{
    (tcrossprod(x, y) * gamma + coef0) ^ degree
}

#' @title linear kernel between X and Y
#' @param x: matrix N rows of samples, P columns of features
#' @param y: matrix M rows of samples, P columns of features
#' @details when Y is is NULL, K(X, X) is computed.
lnr <- function(x, y=NULL, gamma=1/NCOL(x), ...)
{
    tcrossprod(x, y) * gamma
}

cmb <- function(k, w)
{
    if(is.matrix(k))
        k <- list(k)
    w <- as.matrix(w)
    M <- ncol(w)
    L <- nrow(w)
    lapply(1:M, function(m) Reduce(`+`, mapply(`*`, k, w[, m], SIMPLIFY=FALSE)))
}

## Folllowings are the definition of base kernels
ibs <- function(x, level=2, ...)
{
    if(is.null(level))
        x <- apply(x, function(.) . / max(.))
    else
        x <- x / level
    1 - as.matrix(dist(x, 'man')) / ncol(x)
}
