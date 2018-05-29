## Modified MINQUE
Modified_MINQUE_Solver <- function(ViList, X, p)
{
    A        <- LMM_MINQUE_Solver(ViList = ViList, X = X, p = p)$A;
    A        <- (A + t(A)) / 2;
    eigA     <- eigen(A);
    eigA_val <- pmax(eigA$values,0);
    eigA_vec <- eigA$vectors;
    Anew     <- eigA_vec %*% diag(eigA_val) %*% t(eigA_vec);
    
    return(Anew)
}


## MINQUE for LMM
LMM_MINQUE_Solver <- function(ViList, X, p)
{
    n     <- dim(X)[1];
    k     <- length(ViList);
    I     <- diag(rep(1,n));
    S     <- matrix(0, nrow = k, ncol = k);
    AList <- list();
    
    ## Calculate V matrix and its inverse
    V    <- Reduce(f = '+', x = ViList); 
    invV <- chol2inv(chol(V));
    
    ## Calculate P matrix and Q matrix
    Pv   <- X %*% ginv(t(X) %*% invV %*% X) %*% t(X) %*% invV;
    Qv   <- I - Pv;
    
    ## Caculate S matrix to solve for lambda
    for(i in 1:k)
    {
        AList[[i]] <- t(Qv) %*% invV %*% ViList[[i]] %*% invV %*% Qv
    }
    
    for(i in 1:k)
    {
        for(j in 1:k)
        {
            S[i,j] <- sum(AList[[i]] * ViList[[j]]);
        }
    }
    lambda <- ginv(S) %*% p;
    
    ## Calculate A matrix
    A <- mapply(AList, lambda, FUN = '*', SIMPLIFY = FALSE);
    A <- Reduce(f = '+', A);
    
    returnlist <- list(A = A, S = S);
    return(returnlist)
}

knl.mnq <- function(y, K, p=NULL)
{
    ## number of kernels
    k <- length(K)
    ## the contrast matrix to get all variance components
    if(is.null(p))
        p <- diag(k)

    print('begin MINQUE')
    time0 <- proc.time()
    W <- .Call('_knn_knl_mnq_cpp', PACKAGE = 'knn', as.matrix(y), K, p)$f
    time1 <- proc.time()
    print('end MINQUE')

    ## make predictions
    prd <- knl.prd(y, K, W, logged=FALSE)

    ## timing
    rtm <- DF(key='rtm', val=(time1 - time0)['elapsed'])

    ret <- list(par=drop(W), rpt=rbind(rtm, prd))
}    

knl.mnq.R <- function(y, V, P=NULL)
{
    k     <- length(V);
    if(is.null(P))
        P <- diag(k)
    ## I     <- diag(nrow(X))

    print('begin MINQUE')
    time0 <- proc.time()

    ## sum of V_i, i = 1 .. k    # chapter 7
    sumV <- Reduce(f = '+', x = V); 
    
    ## b) Calculate P matrix and Q matrix
    ## Pv = X \(X' \sumV X) X' \sumV   # projection defined by eq 1.2
    ## Qv = I - Pv              # eq 7.1, or lemma 3.4
    ## R = \V Qv                # eq 7.1, for the best choice of A*
    ## R is simplified since X=0 (i.e., not mixed, only kernels)
    R <- chol2inv(chol(sumV));
    
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

    ## lambda(s)
    L <- P %*% MASS::ginv(S)

    ## Calculate A matrix
    A <- list()
    W <- double(nrow(L))
    for(i in 1:nrow(L))
    {
        A[[i]] <- Reduce('+', mapply('*', B, L[i, ], SIMPLIFY = FALSE))
        W[i] <- crossprod(y, A[[i]] %*% y)
    }

    time1 <- proc.time()
    print('end MINQUE')

    ## make predictions
    prd <- knl.prd(y, V, W, logged=FALSE)

    ## timing
    rtm <- DF(key='rtm', val=(time1 - time0)['elapsed'])

    ret <- list(par=drop(W), rpt=rbind(rtm, prd))
    ret
}

pkn.mnq.R <- function(y, V, const = 1, order = 2)
{
    N <- NROW(y)                        # sample size
    vnm <- names(V)                     # kernel names
    if(is.null(vnm))
    {
        vnm <- sprint('v%02d', seq.int(L))
        names(V) <- vnm
    }

    V <- c(list(one=1), V)              # prepend an one
    X <- rep(0, N)                      # no fixed effect
    L <- length(V);                     # number of kernels
    vnm <- names(V)                     # updated names

    ## polynomial kernel expansion
    pke <- gtools::combinations(L, order, repeats.allowed=TRUE)
    K <- list()
    nm <- character(nrow(pke))
    for(k in seq.int(nrow(pke)))
    {
        idx <- pke[k, ]
        K[[k]] <- Reduce(`*`, V[idx])
        nm[k] <- paste(vnm[idx], collapse='.')
    }
    K[[1]] <- diag(N)
    nm[1] <- 'e'
    names(K) <- nm
    k <- length(K);

    ## contrast matrix
    P <- diag(k)

    list(V=V, K=K, y=y, P=P, pke=pke) 
}
