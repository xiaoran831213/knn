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

knl.mnq.R <- function(y, V, P=NULL)
{
    k     <- length(V);
    if(is.null(P))
        P <- diag(k)
    ## I     <- diag(nrow(X))

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

    ## Calculate A matrix and contrasts
    A <- list()
    W <- double(nrow(L))
    for(i in 1:nrow(L))
    {
        A[[i]] <- Reduce('+', mapply('*', B, L[i, ], SIMPLIFY = FALSE))
        W[i] <- crossprod(y, A[[i]] %*% y)
    }

    ## pack up
    list(f=drop(W), A=A)
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
knl.mnq <- function(y, V, order=1, cpp=TRUE)
{
    N <- NROW(y)                        # sample size
    K <- knl.ply(V, order)
    k <- length(K)                      # kernel count
    
    ## call the kernel MINQUE core function
    ## fixed contrast matrix
    P <- diag(k)                        # now k == K.count
    print('begin MINQUE')
    time0 <- proc.time()
    if(cpp)
        W <- .Call('_knn_knl_mnq_cpp', PACKAGE = 'knn', as.matrix(y), K, P)$f
    else
        W <- knl.mnq.R(y, K, P)$f
    time1 <- proc.time()
    print('end MINQUE')

    ## make predictions
    prd <- knl.prd(y, K, W, logged=FALSE, pinv=TRUE)

    ## timing
    rtm <- DF(key='rtm', val=(time1 - time0)['elapsed'])

    ret <- list(par=drop(W), rpt=rbind(rtm, prd))
}

knl.mnq.evl <- function(y, V, vcs, order=1, ...)
{
    K <- knl.ply(V, order)
    knl.prd(y, K, vcs, logged=FALSE, ...)
}
