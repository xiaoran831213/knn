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

knl.mnq <- function(y, K)
{
    L <- length(K)

    ## use MINQUE Solver
    p <- diag(L)

    print('begin MINQUE')
    time0 <- proc.time()
    W <- sapply(seq(L), function(l)
    {
        A <- knl.mnq.A(K, p[, l])$A
        crossprod(y, A %*% y)
    })
    time1 <- proc.time()
    print('end MINQUE')

    ## make predictions
    prd <- knl.prd(y, K, W, logged=FALSE)

    ## timing
    rtm <- DF(key='rtm', val=(time1 - time0)['elapsed'])

    ret <- list(par=W, rpt=rbind(rtm, prd))
}    

knl.mnq.A <- function(v, p)
{
    k     <- length(v);
    ## I     <- diag(nrow(X))
    
    ## sum of V_i, i = 1 .. k    # chapter 7
    V <- Reduce(f = '+', x = v); 
    
    ## b) Calculate P matrix and Q matrix
    ## Pv = X \(X' \V X) X' \V   # projection defined by eq 1.2
    ## Qv = I - Pv              # eq 7.1, or lemma 3.4
    ## R = \V Qv                 # eq 7.1, for the best choice of A*
    ## R is simplified since X=0 (i.e., not mixed, only kernels)
    R <- chol2inv(chol(V));
    
    ## Caculate S matrix, S_ij = Tr(W_i V_j)
    ## B_i = R V_i R
    B <- list()
    for(i in 1:k)
        B[[i]] <- R %*% v[[i]] %*% R
    
    S <- matrix(.0, k, k)
    for(i in 1:k)
    {
        for(j in 2:k)
        {
            S[i, j] <- sum(B[[i]] * v[[j]])
            S[j, i] <- S[i, j]
        }
        S[i, i] <- sum(B[[i]] * v[[i]])
    }
    lambda <- MASS::ginv(S) %*% p
    
    ## Calculate A matrix
    A <- mapply('*', B, lambda, SIMPLIFY = FALSE)
    A <- Reduce('+', A);

    ## return
    list(A=A, S=S, lambda=lambda)
}

## A Helper Function to evaluate r(A,nu) in Infeasible_MINQUE_Newton
evalR <- function(A, nu, Psi, V, t, p, eps)
{
    I        <- diag(rep(1,nrow(A)));
    invA     <- solve(A+eps*I);
    r_dual   <- 2*t*c(V %*% A %*% V) - c(invA) + Psi %*% nu;
    r_primal <- crossprod(Psi, c(A)) - p;
    r        <- c(r_dual, r_primal)
    r_norm   <- sqrt(crossprod(r))
    
    returnList <- list(r_norm = r_norm, r_dual = r_dual, r_primal = r_primal, invA = invA);
    return(returnList)
}

test <- function()
{
    X <- matrix(rnorm(32), 4, 8)
    K <- list(e=diag(nrow(X)), p=tcrossprod(X))
    P <- rbind(diag(length(K)), rep(1, length(K)))
    .Call('_knn_mnq', K, P)
}
