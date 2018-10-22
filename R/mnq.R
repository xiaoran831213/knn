## Modified MINQUE

.mn0 <- function(y, v., w.=NULL, X=NULL, ...)
{
    dot <- list(...)
    psd <- dot$psd %||% 1L              # PSD projection?
    K   <- length(v.)

    ## initial weights
    if(is.null(w.))
        w. <- matrix(1, K, 1)

    ## sum of V_i, i = 1 .. K
    V <- cmb(v., w.)[[1]]
    invV <- try(chol2inv(chol(V)), silent=TRUE)
    if(inherits(invV, 'try-error'))
        invV <- solve(V)
    
    ## b) Calculate P matrix and Q matrix
    ## Pv = X \(X' V^{-1} X) X' V^{-1}   # projection defined by eq 1.2
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
    
    ## B_i = R V_i R
    B <- list()
    for(i in 1:K)
        B[[i]] <- R %*% v.[[i]] %*% R

    ## Caculate F: F_ij = Tr(R V_i R V_j)
    F <- matrix(.0, K, K)
    for(i in 1:K)
    {
        for(j in 2:K)
        {
            F[i, j] <- sum(B[[i]] * v.[[j]])
            F[j, i] <- F[i, j]
        }
        F[i, i] <- sum(B[[i]] * v.[[i]])
    }

    ## lambda_i = P_i * F^{-1}, to solve for each variance compoent, 
    ## the contract matrix P must be I(k), thus Lambda == F^{-1}.
    ## solve F L = [p1, p2, ...]
    ## when the K components are desired, [p1, p2, ...]  = I_K.
    ## ==> L = solve(F, I_K) = F^{-}
    L <- MASS::ginv(F)

    ## Calculate A matrix and contrasts
    A <- list()
    w. <- double(nrow(L))               # variance component estimate
    for(i in 1:nrow(L))
    {
        a <- Reduce('+', mapply('*', B, L[i, ], SIMPLIFY = FALSE))
        w <- sum(y * a %*% y)          # attempt to estimate s^2
        if(w < 0 && psd)               # make A matrix PSD if s^2 < 0.
        {
            ## Modified y'A y
            a <- eigen((a + t(a)) / 2)
            d <- pmax(a$values, 0)
            v <- a$vectors
            a <- v %*% (d * t(v))       # A_hat_i: modified A matrix
            w <- sum(y * a %*% y)
        }
        A[[i]] <- a
        w.[i] <- w
    }

    ## estimate standard error of variance components estimate
    C <- cmb(v., w.)[[1]]                 # marginal covariance of y
    e <- double(nrow(L))
    ## var(vc_i) = var(Y' \hat{A}_i y) = 2Tr(\hat{A}_i C \hat{A}_i C)
    for(i in 1:nrow(L))                 # 
        e[i] <- 2 * sum((A[[i]] %*% C)^2)

    ## pack up
    list(vcs=w., se2=e)
}

.ginv <- function(x, e=sqrt(.Machine$double.eps))
{
    with(svd(x),
    {
        i <- d > d[1L] * e
        if(all(i))
            v %*% (t(u) / d)
        else
            v[, i, drop=FALSE] %*% (t(u[, i, drop=FALSE]) / d[i])
    })
}
    
.mnq <- function(y, v., w.=NULL, X=NULL, ...)
{
    dot <- list(...)
    K   <- length(v.)

    ## initial weights
    w. <- w. %||% matrix(1, K, 1)

    ## sum of V_i, i = 1 .. K
    V <- cmb(v., w., TRUE)
    
    ## b) Calculate Q, R V_i, and R y, avoid calculating R directly.
    if(!is.null(X))
    {
        . <- solve(V, X)     # V^{-1}X
        P <- X %*% ginv(t(X) %*% .) %*% .

        ## R V_i = V^{-1}(I - P) V_i = V^{-1}(V_i - P V_i)
        Rv. <- list()
        for(i in seq.int(K))
            Rv.[[i]] <- solve(V, v.[[i]] - P %*% v.[[i]])

        ## R y   = V^{-1}(I - P) y   = V^{-1}(y   - P y  )
        Ry <- solve(V, y - P %*% y)
    }
    else
    {
        ## R V_i = V^{-1} V_i
        Rv. <- list()
        for(i in seq.int(K))
        {
            print(i)
            Rv.[[i]] <- solve(V, v.[[i]])
        }
        ## R y   = V^{-1} y
        Ry <- solve(V, y)
    }
    ## v_i = e' V_i e = y' R V_i R y
    v_<- double(K)
    for(i in seq.int(K))
        v_[i] <- sum(Ry * v.[[i]] %*% Ry)

    ## Caculate F: F_ij = Tr(R V_i R V_j)
    F <- matrix(.0, K, K)
    for(i in seq.int(K))
    {
        for(j in seq.int(2L, K))
        {
            F[i, j] <- sum(Rv.[[i]] * Rv.[[j]])
            F[j, i] <- F[i, j]
        }
        F[i, i] <- sum(Rv.[[i]] * Rv.[[i]])
    }

    ## [s2_1, s2_2, ..., s2_K]^T =
    ## [yA1y, yA2y, ..., yAKy]^T = solve(F, v_)
    w <- solve(F, v_)
    ## [ A1, A2, ... AK] is not directly calculated

    ## pack up
    list(vcs=w, se2=NULL)
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
knl.mnq <- function(y, V, W=NULL, cpp=TRUE, itr=1, ...)
{
    dot <- list(...)
    prd <- dot$prd %||% FALSE           # make self prediction
    psd <- dot$psd %||% 1L              # make PSD projection?
    zbd <- dot$zbd %||% 1L              # 0 as lower bound?
    tol <- dot$tol %||% 1e-5
    N <- NROW(y)                        # sample size
    Q <- NCOL(y)                        # response count
    hst <- dot$hst %||% list()          # history

    ## append noise kernel
    C <- c(list(eps=diag(N)), V)
    K <- length(C)                      # kernel count
    if(is.null(W))                      # initial weights
        W <- matrix(1, K, Q)

    ## print('begin MINQUE')
    t0 <- Sys.time()

    while(itr > 0)
    {
        ## call the kernel MINQUE core function
        if(cpp)
            r <- .Call('knl_mnq', as.matrix(y), C, psd, PACKAGE='knn')
        else
            r <- .mnq(y, C, W, psd=psd)
        Z <- r$vcs
        if(zbd)                         # zero bounded?
            Z <- pmax(Z, 0)

        D <- max(abs(Z - W))
        cat('MNQ:', sprintf('%03d', itr), c(Z, D), '\n')
        if(D < tol)
            break
        W <- Z
        itr <- itr - 1
    }
    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    ## print('end MINQUE')

    rpt <- rbind(vpd(y, V, Z), rtm=DF(key='rtm', val=td))
    names(Z) <- names(C)

    list(par=Z, se2=r$se2, rpt=rpt)
}
