## Modified MINQUE

knl.psd <- function(A)
{
    X <- (A + t(A)) / 2
    with(eigen(X), list(v=vectors, d=pmax(values, 0)))
}

knl.mnq.R <- function(y, V, X=NULL)
{
    k     <- length(V)

    ## sum of V_i, i = 1 .. k    # chapter 7
    sumV <- Reduce(f = '+', x = V)
    invV <- chol2inv(chol(sumV))
    
    ## b) Calculate P matrix and Q matrix
    ## Pv = X \(X' \sumV X) X' \sumV   # projection defined by eq 1.2
    ## Qv = I - Pv              # eq 7.1, or lemma 3.4
    ## R = \V Qv                # eq 7.1, for the best choice of A*
    ## R is simplified since X=0 (i.e., not mixed, only kernels)
    if(!is.null(X))
    {
        I <- diag(nrow(X))
        . <- crossprod(X, invV)
        R <- invV %*% (I - X %*% MASS::ginv(. %*% X) %*% .)
    }
    else
        R <- invV
    
    ## Caculate S matrix, S_ij = Tr(W_i V_j)
    ## B_i = R V_i R
    B <- list()
    for(i in 1:k)
        B[[i]] <- R %*% V[[i]] %*% R

    S <- matrix(.0, k, k)
    for(i in 1:k)
    {
        for(j in 2:k)
        {
            S[i, j] <- sum(B[[i]] * V[[j]])
            S[j, i] <- S[i, j]
        }
        S[i, i] <- sum(B[[i]] * V[[i]])
    }

    ## lambda_i = P_i * S^{-1}, to solve for each variance compoent, 
    ## the contract matrix P must be I(k), thus Lambda == S^{-1}.
    L <- MASS::ginv(S)

    ## Calculate A matrix and contrasts
    A <- list()
    W <- double(nrow(L))                # variance component estimate
    for(i in 1:nrow(L))
    {
        a <- Reduce('+', mapply('*', B, L[i, ], SIMPLIFY = FALSE))
        w <- sum(y * a %*% y)          # attempt to estimate s^2
        if(w < 0)                      # make A matrix PSD if s^2 < 0.
        {
            ## Modified y'A y
            a <- eigen((a + t(a)) / 2)
            d <- pmax(a$values, 0)
            v <- a$vectors
            a <- v %*% (d * t(v))       # A_hat_i: modified A matrix
            w <- sum(y * a %*% y)
        }
        A[[i]] <- a
        W[i] <- w
    }

    ## estimate standard error of variance components estimate
    C <- cmb(V, W)[[1]]                 # marginal covariance of y
    e <- double(nrow(L))
    ## var(vc_i) = var(Y' \hat{A}_i y) = 2Tr(\hat{A}_i C \hat{A}_i C)
    for(i in 1:nrow(L))                 # 
        e[i] <- 2 * sum((A[[i]] %*% C)^2)

    ## pack up
    list(vcs=W, se2=e)
}

#' Kernel Polynomial Expansion
#'
#' @param V list of kernels (matrices of covariance)
#' @param order r the order of polynomial
#' @param v index of the basic kernels
#'
#' @return the polynomial terms formed by Schur products of the basic kernels,
#' up to the specified order, including a zero order identidy kernel for noise.
knl.ply <- function(V, order=1)
{
    L <- length(V)                      # kernel count
    N <- names(V)                       # kernel names
    if(is.null(N))
        N <- letters[seq.int(L)]

    V <- c(list(I=1), V)                # prepend "1"
    N <- c('.', N)                      # updated names
    L <- length(V)                      # updated count

    .kp <- function(n, r, v = 1:n)      # see gtools::combinations
    {
        v0 <- vector(mode(v), 0)
        if (r == 0) 
            v0
        else if (r == 1) 
            matrix(v, n, 1)
        else if (n == 1) 
            matrix(v, 1, r)
        else
            rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 1, r, v[-1]))
    }

    P <- .kp(L, order)                  # polynomial indicies
    K <- list()                         # expended kernels
    M <- character(nrow(P))             # expanded kernel names
    for(k in 1:nrow(P))
    {
        K[[k]] <- Reduce(`*`, V[P[k, ]])
        M[k] <- paste(N[P[k, ]], collapse='')
    }

    K[[1]] <- diag(NROW(K[[2]]))       # replace "1" with identity kernel
    names(K) <- M                      # assign names
    K
}

#' Kernel MINQUE
#'
#' @param y a vector of dependent variable
#' @param V a list of kernels (matrices of covariance), with no identity
#' kernel (the noise term).
#' @param order an interger controls the kernel polynomial expansion
#' @param cpp use compiled Cpp code instead of R routine.
#' @return a list of statistics:
#' 
#' * par: a vector variance component esimates
#' * a table of performance measures
#'   - running time;
#'   - mean square error between y and y.hat;
#'   - negative log likelihood assuming y ~ N(0, sum_i(V_i * par_i))
#'   - cor(y, y.hat)
knl.mnq <- function(y, V, order=1, cpp=TRUE, ...)
{
    dot <- list(...)
    prd <- dot$prd %||% FALSE           # make self prediction
    psd <- dot$psd %||% 1L              # make PSD projection?
    
    N <- NROW(y)                        # sample size

    ## print('begin MINQUE')
    t0 <- Sys.time()
    C <- knl.ply(V, order)
    k <- length(C)                      # kernel count
    
    ## call the kernel MINQUE core function
    ## fixed contrast matrix
    if(cpp)
        r <- .Call('knl_mnq', PACKAGE = 'knn', as.matrix(y), C, psd)
    else
        r <- knl.mnq.R(y, C)
    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    ## print('end MINQUE')

    ## make predictions
    if(prd)
        prd <- vpd(y, V, r$vcs, ln=0, rt=1)
    else
        prd <- NULL

    ## timing
    rtm <- DF(key='rtm', val=td)
    ret <- list(par=drop(r$vcs), se2=drop(r$se2), rpt=rbind(rtm=rtm, prd))
}
