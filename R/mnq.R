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

pkn.mnq.R <- function(y, V, P, const = 1, order = 2)
{
    L <- length(V);                     # number of kernels
    N <- NROW(y)
    IdMat <- diag(N);                   # identity kernel
    X <- rep(0, N)                      # no fixed effect

    terms <- do.call(expand.grid, replicate(3, 0L:1L, simplify=FALSE))
    ## powerMat <- perm(order, L+1);
    ## coef     <- powerMat[,ncol(powerMat)];
    ## powerMat <- powerMat[,1:(ncol(powerMat)-1)];
    VList    <- list();
    
    newBaseKernelList <- append(V, const, 0);
    
    for(k in 1:nrow(powerMat))
    {
        tempList   <- mapply(newBaseKernelList, FUN = '^', powerMat[k,], SIMPLIFY = FALSE);
        VList[[k]] <- Reduce(f = '*', tempList);
        ## VList      <- mapply(VList, FUN = '*', coef, SIMPLIFY = FALSE);
    }
    
    VList[[k + 1]] <- IdMat;
    
    k <- length(VList);
    minque_est <- rep(0,k);
    
    for(i in 1:k)
    {
        P <- rep(0,k);
        P[i] <- 1;
        minque_mat <- LMM_MINQUE_Solver(ViList = VList, X = X, p = P);
        minque_est[i] <- t(y) %*% minque_mat %*% y;
    }

    returnList <- list(minque_est = minque_est, VList = VList)
    return(returnList)
}
