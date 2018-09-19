.mq <- function(y, knl, W=NULL, X=NULL, psd=FALSE)
{
    K <- length(knl)
    N <- NROW(y)
    Q <- NCOL(y)
    if(is.null(W))
        W <- matrix(1/K, K, Q)

    ## sum of V_i, i = 1 .. k    # chapter 7
    V <- cmb(knl, W)[[1]]
    Z <- chol2inv(chol(V))
    
    ## b) Calculate P matrix and Q matrix
    ## Pv = X [(X'V^{-1} X)^{-1}] X'V^{-1}   # projection defined by eq 1.2
    ## Qv = I - Pv              # eq 7.1, or lemma 3.4
    ## R = \V Qv                # eq 7.1, for the best choice of A*
    ## R is simplified if X=0 (i.e., no fix effect, only random)
    if(!is.null(X))
    {
        Pv <- X %*% ginv(t(X) %*% Z %*% X) %*% t(X) %*% Z
        ## solve(t(X) %*% solve(V, X), t(X))
        Qv <- diag(N) - Pv
        R <- Z %*% P
    }
    else
        R <- Z
    
    ## Caculate S matrix, S_ij = Tr(W_i V_j)
    ## B_i = R V_i R
    B <- list()
    for(i in 1:K)
        B[[i]] <- R %*% knl[[i]] %*% R

    S <- matrix(.0, K, K)
    for(i in 1:K)
    {
        for(j in 2:K)
        {
            S[i, j] <- sum(B[[i]] * knl[[j]])
            S[j, i] <- S[i, j]
        }
        S[i, i] <- sum(B[[i]] * knl[[i]])
    }

    ## lambda_i = P_i * S^{-1}, to solve for each variance compoent, 
    ## the contrast matrix P must be I(k), thus Lambda == S^{-1}.
    Lambda <- ginv(S)

    ## Calculate A matrix and contrasts
    W <- double(K)                # variance component estimate
    for(i in 1L:K)
    {
        a <- Reduce('+', mapply('*', B, Lambda[i, ], SIMPLIFY=FALSE))
        w <- sum(y * a %*% y)
        if(w < 0 && psd)                      # make A matrix PSD if s^2 < 0.
        {
            ## Modified y'A y
            a <- eigen((a + t(a)) / 2)
            d <- pmax(a$values, 0)
            v <- a$vectors
            a <- v %*% (d * t(v))       # A_hat_i: modified A matrix
            w <- sum(y * a %*% y)
        }
        W[i] <- w
    }

    ## no attempt to estimate s^2, since A[[i]] may be invalid
    ## pack up
    list(vcs=W, se2=NA)
}

itr.mnq <- function(y, knl, W=NULL, X=NULL, psd=TRUE)
{
    w0 <- .mq(y, knl, NULL, X, psd=FALSE)$vcs
    pa <- list()
    for(i in 1:50)
    {
        print(w0)
        w1 <- .mq(y, knl, w0, X, psd=psd)$vcs
        if(max(abs(w0 - w1)) < 1e-5)
            break
        w0 <- w1
        pa <- c(pa, list(w0))
    }
    w0 <- mean(pa)
    list(par=w0, se2=NA)
}
