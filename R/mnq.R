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
    P <- cbind(diag(length(K)), rep(1, length(K)))
    .Call('_cpp_mnq', K, P)
}
