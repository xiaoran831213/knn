## Modified MINQUE

.mn0 <- function(y, v., w.=NULL, X=NULL, ...)
{
    dot <- list(...)
    psd <- dot$psd %||% 0L              # PSD projection?
    K   <- length(v.)

    ## initial weights
    if(is.null(w.))
        w. <- matrix(1, K, 1)

    ## sum of V_i, i = 1 .. K
    V <- cmb(v., w., TRUE)
    ## invV <- try(chol2inv(chol(V)), silent=TRUE)
    ## if(inherits(invV, 'try-error'))
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
        R <- invV %*% (I - X %*% .ginv(. %*% X) %*% .)
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
    L <- .ginv(F)

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
    ## C <- cmb(v., w., TRUE)              # marginal covariance of y
    ## e <- double(nrow(L))
    ## var(vc_i) = var(Y' \hat{A}_i y) = 2Tr(\hat{A}_i C \hat{A}_i C)
    ## for(i in 1:nrow(L))                 # 
    ##     e[i] <- 2 * sum((A[[i]] %*% C)^2)

    ## pack up
    names(w.) <- names(v.)
    list(vcs=w., se2=NULL)
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
    w. <- w. %||% rep(1, K)

    ## sum of V_i, i = 1 .. K
    V <- cmb(v., w., TRUE)
    
    ## Calculate P, and Q = I - P, then calculate R V_i and R y
    ## avoid direct calculation of R = V^{-1} Q
    Rv. <- list()
    if(is.null(X))
    {
        ## R V_i = V^{-1}(I - P) V_i = V^{-1}(V_i - P V_i)
        for(i in seq.int(K))
            Rv.[[i]] <- solve(V, v.[[i]])

        ## R y   = V^{-1}(I - P) y   = V^{-1}(y   - P y  )
        Ry <- solve(V, y)
    }
    if(!is.null(X))
    {
        B <- solve(V, X)                # V^{-1}X
        P <- X %*% .ginv(t(X) %*% B) %*% t(B)

        ## R V_i = V^{-1}(I - P) V_i = V^{-1}(V_i - P V_i)
        for(i in seq.int(K))
            Rv.[[i]] <- solve(V, v.[[i]] - P %*% v.[[i]])

        ## R y   = V^{-1}(I - P) y   = V^{-1}(y   - P y  )
        Ry <- solve(V, y - P %*% y)
    }

    ## v_i = e' V_i e = y' R V_i R y
    v_<- double(K)
    for(i in seq.int(K))
        v_[i] <- sum(Ry * v.[[i]] %*% Ry)

    ## Caculate F: F_ij = Tr(R V_i R V_j)
    F <- matrix(.0, K, K)
    for(i in seq.int(K))
    {
        for(j in seq.int(K)[-1L])
        {
            F[i, j] <- sum(Rv.[[i]] * Rv.[[j]])
            F[j, i] <- F[i, j]
        }
        F[i, i] <- sum(Rv.[[i]] * Rv.[[i]])
    }

    ## [s2_1, s2_2, ..., s2_K]^T =
    ## [yA1y, yA2y, ..., yAKy]^T = solve(F, v_) = pinv(F) v_
    ## w <- solve(F, v_)
    w <- .ginv(F) %*% v_ 
    ## [ A1, A2, ... AK] is not directly calculated

    ## GLS for fixed effects
    if(!is.null(X))
    {
        b <- .ginv(crossprod(X, B)) %*% crossprod(B, y)
        names(b) <- colnames(X)
    }
    else
        b <- NULL

    ## pack up
    list(vcs=w, se2=NULL, fix=b)
}

#' Generalized Least Square estimates of beta
#'
#' \hat{beta} = (X^T V^{-1} X)^{+} X^T V^{-1} y
#' 
#' @param V matrix of covariance
#' @param X matrix of covariate
#' @param y vector of response
.gls <- function(V, X, y)
{
    ## V^{-1} X
    Z <- solve(V, X)

    ## (X^T Z)^{+} Z^T y
    ## .ginv(t(X) %*% Z) %*% t(Z) %*% y
    b <- .ginv(crossprod(X, Z)) %*% crossprod(Z, y)
    names(b) <- colnames(X)
    b
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
#' @param itr number of iterations, def=10
#' @param fix estimate fixed effect, fef=0
#' @return a list of statistics:
#' 
#' * par: a vector variance component esimates
#' * a table of performance measures
#'   - running time;
#'   - mean square error between y and y.hat;
#'   - negative log likelihood assuming y ~ N(0, sum_i(V_i * par_i))
#'   - cor(y, y.hat)
knl.mnq <- function(y, V, X=NULL, W=NULL, cpp=FALSE, itr=30, ...)
{
    dot <- list(...)
    psd <- dot$psd %||% 0L              # make PSD projection?
    zbd <- dot$zbd %||% 0L              # zero as lower bound?
    vbd <- dot$vbd %||% 1L              # upper bound considered?
    ebd <- dot$ebd %||% 1L              # bound Epsilon always?
    tol <- dot$tol %||% 1e-5
    vth <- dot$vth %||% 0.1             # variance threshold

    ## append noise kernel
    N <- length(y)                      # sample size
    C <- c(list(EPS=diag(N)), V)
    K <- length(C)                      # kernel count
    nm0 <- names(C)

    ## prepend intercept if necessary
    if(is.null(X) || !any(grepl('X00', names(X)[1], TRUE)))
        X <- cbind(X00=rep(1.0, N), X)

    ## default total variance without fixed effect
    SST <- sum(residuals(lm(y ~ X))^2) / (N-1)  # sum(y^2) / (N-1)

    ## initial weights
    ini <- if(is.null(W)) rep(SST/K, K) else W
    vcs <- ini
    par <- ini

    ## kappa of kernels
    kpa <- c(1, sapply(V, kappa, method="direct")) # initial kappa
    msk <- c(EPS=Inf, rank(kpa[-1]))               # kernel mask

    ## print('begin MINQUE')
    t0 <- Sys.time()

    ps0 <- function(...) paste(..., collapse=" ")
    sp0 <- function(...) ps0(sprintf(...))
    ca0 <- function(...) cat(..., "\n", sep="")
    hdr <- ps0('MN:', sp0('%3s', "ITR"), sp0("%6s", nm0), sp0("%7s", "MISC"))

    ## helper inner functions
    .df <- function() max(abs(vcs - par), na.rm=TRUE)
    .DF <- function() sp0("%6.4fD", .df())

    .ev <- function()
    {
        v <- c(SUM=sum(vcs, na.rm=TRUE), vcs[!is.na(vcs)])
        max(ifelse(v < 0, -v / SST, v / SST - 1))
    }
    .EV <- function() sp0("%6.4fV", .ev())
    
    VC <- function() ps0('VC:', sp0('%03d', itr), sp0("%6.3f", vcs), .DF(), .EV())

    ## print the header
    ca0(hdr)

    ca0(VC())
    while(itr > 0)
    {
        par <- vcs
        ## call MINQUE core function, skip kernels with zero weights.
        if(cpp)
            r <- .Call('knl_mnq', as.matrix(y), C[!is.na(vcs)], psd, PACKAGE='knn')
        else
            r <- .mnq(y, C[!is.na(vcs)], par[!is.na(vcs)], X=X)

        ## total variance (conditioned on fixed effect)
        if(!is.null(r$fix))
            SST <- sum((y - X %*% r$fix)^2) / (N-1)

        ## get solution, and zero bounded?
        vcs[!is.na(vcs)] <- r$vcs
        if(zbd)
            vcs[!is.na(vcs) & vcs < 0] <- 0.0

        ## check excessive error
        if(ebd && vcs[1] < 0)
        {
            msg <- ps0('VC:', sp0('%03d', itr), sp0("%6.3f", vcs), "   EBD")
            ca0(msg)
            msk  <- msk - 1
            vcs[msk <= 0] <- NA
            vcs[!is.na(vcs)] <- ini[!is.na(vcs)]
            next
        }

        ## check excessive variance component
        if(vbd)
        {
            tmp <- na.omit(vcs)
            tmp <- c(SUM=sum(tmp), tmp) / SST
            tmp <- ifelse(tmp < 0, abs(tmp), tmp - 1)
            ## print(max(tmp))
            tmp <- max(tmp)
            if(tmp > vth && sum(!is.na(vcs)) > 1)
            {
                msg <- ps0('VC:', sp0('%03d', itr), sp0("%6.3f", vcs), "   VBD", sp0("%6.4fV", tmp))
                ca0(msg)
                msk <- msk - 1
                vcs[msk <= 0] <- NA
                vcs[!is.na(vcs)] <- ini[!is.na(vcs)]
                next
            }
        }

        ## convergence check
        if(sum(is.na(vcs)) == sum(is.na(par)))
            dff <- max(abs(vcs - par), na.rm=TRUE)
        else
            dff <- Inf
        ca0(VC())
        if(dff < tol)
            break
        itr <- itr - 1
    }
    vcs[is.na(vcs)] <- 0

    ## fixed effects
    if(!is.null(X))
        B <-.gls(cmb(C, vcs, TRUE), X, y)
    else
        B <- NULL
    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    ## print('end MINQUE')

    ## pack up and return: estimates
    names(vcs) <- nm0
    par <- c(B, vcs)
    rpt <- c(vpd(par, y, V, X, rt=0))

    ## reports & return
    ret <- list(par=par, kpa=kpa, rpt=rpt, rtm=td)
    ret
}
