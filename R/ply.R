## polynomial kernels

#' polynomial kernel
#'
#' @param ... a series of basic kernel functions. A kernel function can
#' either be a function, or a list pre-calculated kernel matrices.
#' @param D   integer, degree of the polynomial expansion
#' @param O   integer, orthogonality
#' 
#' * 0: ordinary polynomial
#' * 1: orthogonal
#' * 2: priciple component
#' 
#' @param J: integer, use non-zero to include a J matrix of all ones.
#' 
#' @return a function to eveluate a series of N by N kernel matrices
#' from an input of N samples and P features.
PK <- function(..., D=2, O=0, J=0)
{
    ev <- environment()
    cl <- match.call()
    ## make a functional that builds high order kernel(s) from x
    structure(function(x)
    {
        ## either apply the basic functions to x, or keep the
        ## pre-calculated kernel matrices as it is.
        B <- list(...)
        B <- lapply(B, function(b)
        {
            if(is.function(b))
                b <- b(x)
            b
        })
        B <- unlist(B, FALSE)
        A <- B[[1]]                     # a template matrix
        utr <- upper.tri(A, TRUE)
        B <- lapply(B, `[`, utr)
        nms <- names(B)
        
        ## expansion
        arg <- c(B, list(degree=D, coefs=ev[['oc']], raw= (O != 1)))
        B <- do.call(polym, arg)
        
        ## keep orthgonal coeficients
        oc <- attr(B, 'coefs')
        ev[['oc']] <- if(length(nms) >1) oc else if(is.null(oc)) NULL else list(oc)
        ## print(ev[['oc']])

        ## pca?
        if(O == 2)
        {
            pc <- ev[['pc']]
            if(is.null(pc))
            {
                print("New PCA.")
                pc <- prcomp(B, retx=FALSE, scale.=TRUE)
                ev[['pc']] <- pc
            }
            else
                print(pc)
            B <- predict(pc, B)
            nms <- NULL
        }
        
        B <- as.data.frame(B)
        names(B) <- lapply(strsplit(names(B), '[.]'), function(pls)
        {
            n <- nms[pls != "0"]
            p <- pls[pls != "0"]
            paste0(n, p, collapse='.')
        })

        ## restore the kernel matrices
        B <- lapply(B, function(k)
        {
            A[utr] <- k; A <- t(A); A[utr] <- k; A
        })
        if(J)
        {
            B <- c(JX1=unname(JX(x)), B)
        }
        B
    }, class=c('pk', 'function'))
}

#' Print a Polynomial Kernel Object
print.pk <- function(x, ...) print(get('cl', environment(x)), ...)

#' Run Simulation for a Polynomial Object
simulate.pk <- function(x, N=5, P=8, ...)
{
    u <- matrix(rnorm(N * P), N, P)
    lapply(x(u), round, 4)
}

test.pk <- function()
{
    N=5
    P=8
    MX <- list(
        A1=tcrossprod(matrix(rnorm(N * P), N, P)),
        A2=tcrossprod(matrix(rnorm(N * P), N, P)))

    f <- PK(LN, I2, MX, J=1)
    print(f)
    simulate(f, N, P)
}
