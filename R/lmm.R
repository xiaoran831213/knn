## LMM or pure variance component
lmm.dv1 <- function(W, K, Y, ...)
{
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- exp(W)
    
    V <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(V, chol)                # decompotion
    A <- lapply(H, chol2inv)            # acc: a_i, (i = 1 .. M)

    Y <- as.matrix(Y)

    ## alpha_i = a_i %*% y_i
    YAY <- vapply(1:M, function(m) A[[m]] %*% Y[, m], Y[, 1])

    ## container of first derivative
    dv1 <- matrix(.0, L, M)

    ## tmp_i = alpha_i alpha_i' - a_i
    TMP <- lapply(1:M, function(m) tcrossprod(YAY[, m]) - A[[m]])
    for(m in seq.int(M))
    {
        for(l in seq.int(L))
        {
            ## dvt(y_m, w_l) = -1/2 Tr[(alpha_i alpha_i' - a_i) K_l]
            dv1[l, m] <- -.5 * sum(TMP[[m]] * K[[l]])
        }
    }
    dv1 <- dv1 * W
    dv1
}

lmm.dv2 <- function(W, K, Y, ...)
{
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- exp(W)

    V <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(V, chol)                # decompotion
    A <- lapply(H, chol2inv)            # acc: a_i, (i = 1 .. M)

    Y <- as.matrix(Y)

    ## alpha_i = (a_i y_i) (a_i y_i)'
    YAY <- lapply(1:M, function(m) tcrossprod(A[[m]] %*% Y[, m]))

    ## container of second derivative
    dv2 <- lapply(1:M, function(m)
    {
        dmx <- matrix(.0, L, L)
        for(r in seq.int(1, L))
        {
            for(s in seq.int(r, L))
            {
                tmp <- 2 * YAY[[m]] - A[[m]]
                dmx[r, s] <- .5 * W[r, m] * sum(A[[m]] %*% K[[r]] * tmp %*% K[[s]]) * W[s, m]
                ## dmx[r, s] <- .5 * W[r, m] * tr.hut(A[[m]], K[[r]], tmp, K[[s]], N=1000) * W[s, m]
                dmx[s, r] <- dmx[r, s]
            }
            tmp <- YAY[[m]] - A[[m]]
            dmx[r, r] <- dmx[r, r] -.5 * sum(tmp * K[[r]]) * W[r, m]
        }
        dmx
    })
    dv2
}

lmm.fsi <- function(W, K, Y, ...)
{
    Y <- as.matrix(Y)
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- exp(W)

    V <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(V, chol)                # decompotion
    A <- lapply(H, chol2inv)            # acc: a_i, (i = 1 .. M)

    ## container of second derivative
    dv2 <- lapply(1:M, function(m)
    {
        dmx <- matrix(.0, L, L)
        for(r in seq.int(1, L))
        {
            for(s in seq.int(r, L))
            {
                dmx[r, s] <- .5 * W[r, m] * sum(A[[m]] %*% K[[r]] * A[[m]] %*% K[[s]]) * W[s, m]
                dmx[s, r] <- dmx[r, s]
            }
        }
        dmx
    })
    dv2
}

knl.prd <- function(y, K, W, logged=TRUE, pinv=FALSE)
{
    ## is W logged
    if(logged)
        W <- exp(W)

    ## make predictions
    v <- cmb(K, W)[[1]]
    f <- v - diag(W[1], length(y))      # W[1] is PHI

    ## prediction 1: conditional Gaussian
    h <- f %*% (solve(v) %*% y)
    
    mse <- mean((y - h)^2)
    cyh <- cor(y, h)
    nlk <- nlk(y, v)
    rpt <- DF(key=c('mse', 'nlk', 'cyh'), val=c(mse, nlk, cyh))

    rpt
}

rop.lmm <- function(y, K, W=NULL)
{
    L <- length(K)
    if(is.null(W))
        W <- matrix(rnorm(L, sd=.05), L, 1L)

    obj <- function(x) nlk(y, cmb(K, exp(x))[[1]])
    grd <- function(x) lmm.dv1(x, K, y)
    hsn <- function(x) lmm.dv2(x, K, y)
    fsi <- function(x) lmm.fsi(x, K, y)

    ## using R's optimizer
    ## library(numDeriv)
    ## print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))
    dv1 <- lmm.dv1(W, K, y)             # gradient before
    PF("DV1.MAX = %9.3f\n", dv1[which.max(abs(dv1))])
    time0 <- proc.time()
    opt <- optim(W, obj, grd, method="L-BFGS-B", control=list(trace=1))
    time1 <- proc.time()
    W <- opt$par
    ## print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))
    dv1 <- lmm.dv1(W, K, y)             # gradient after
    PF("DV1.MAX = %9.3f\n", dv1[which.max(abs(dv1))])

    ## make predictions
    prd <- knl.prd(y, K, W)

    ## timing
    rtm <- DF(key='rtm', val=(time1 - time0)['elapsed'])

    opt$rpt <- rbind(rtm, prd)
    opt
}

mnq.lmm <- function(y, K)
{
    L <- length(K)

    ## use MINQUE Solver
    X <- as.matrix(rep(0, NROW(y)))
    p <- diag(L)
    print('begin MINQUE')
    time0 <- proc.time()
    W <- sapply(seq(L), function(l)
    {
        A <- LMM_MINQUE_Solver(K, X, p[, l])$A
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

lmm <- function(y, v, e)
{
    N <- nrow(v)
    f <- v - diag(e, N)
    u <- chol(v)
    a <- chol2inv(u)

    ## prediction 1: conditional Gaussian
    h <- f %*% (a %*% y)
    mse <- mean((y - h)^2)
    cyh <- cor(y, h)

    ## prediction 2: leave one out CV
    ## h <- y - a %*% y / diag(a)
    ## loo <- mean((y - h)^2)

    ## negative log likelihood
    nlk <- nlk(y, v, u, a)

    DF(key=c('mse', 'nlk', 'cyh'), val=c(mse, nlk, cyh))
}
