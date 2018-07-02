## LMM or pure variance component
lmm.dv1 <- function(W, K, Y, ...)
{
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- exp(W)
    
    C <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(C, chol)                # decompotion
    A <- lapply(H, chol2inv)            # acc: a_i, (i = 1 .. M)

    Y <- as.matrix(Y)

    ## alpha_i = a_i %*% y_i
    YAY <- lapply(1:M, function(m) A[[m]] %*% Y[, m])

    ## tmp_i = alpha_i alpha_i' - a_i
    TMP <- lapply(1:M, function(m) tcrossprod(YAY[[m]]) - A[[m]])

    ## 1st derivative
    dv1 <- lapply(seq.int(M), function(m)
    {
        sapply(seq.int(L), function(l)
        {
            -.5 * sum(TMP[[m]] * K[[l]])
        }) * W
    })
    dv1
}

lmm.dv2 <- function(W, K, Y, ...)
{
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- exp(W)

    C <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(C, chol)                # decompotion
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

    C <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(C, chol)                # decompotion
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

knl.prd <- function(y, K, W, logged=TRUE, ...)
{
    ## is W logged
    if(logged)
        W <- exp(W)

    ## make predictions
    v <- cmb(K, W)[[1]]
    f <- v - diag(W[1], length(y))      # W[1] is PHI

    ## prediction 1: conditional Gaussian
    y <- unname(drop(as.vector(y)))
    h <- unname(drop(f %*% (solve(v) %*% y)))

    mse <- mean((y - h)^2)
    ## gurad against zero-standard deviation
    cyh <- tryCatch(cor(y, h), warning=function(w) 0, error=function(e) NA)
    nlk <- nlk(y, v)
    rpt <- DF(key=c('mse', 'nlk', 'cyh'), val=c(mse, nlk, cyh))

    rpt
}

rop.lmm <- function(y, K, W=NULL)
{
    L <- length(K)
    if(is.null(W))
        W <- matrix(rnorm(L, sd=.5), L, NCOL(y))

    obj <- function(x) nlk(y, cmb(K, exp(x))[[1]])
    grd <- function(x) lmm.dv1(x, K, y)[[1]]
    hsn <- function(x) lmm.dv2(x, K, y)[[1]]
    fsi <- function(x) lmm.fsi(x, K, y)

    ## using R's optimizer
    library(numDeriv)
    ## print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))
    dv1 <- lmm.dv1(W, K, y)[[1]]        # gradient before
    PF("DV1.MAX = %9.3f\n", dv1[which.max(abs(dv1))])
    time0 <- proc.time()
    opt <- optim(W, obj, grd, method="L-BFGS-B", control=list(trace=1))
    time1 <- proc.time()
    W <- opt$par
    print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))
    dv1 <- lmm.dv1(W, K, y)[[1]]        # gradient after
    PF("DV1.MAX = %9.3f\n", dv1[which.max(abs(dv1))])

    ## make predictions
    prd <- knl.prd(y, K, W)

    ## timing
    rtm <- DF(key='rtm', val=unname((time1 - time0)['elapsed']))

    opt$rpt <- rbind(rtm, prd)
    opt
}

nwt.lmm <- function(y, K, W=NULL)
{
    L <- length(K)
    if(is.null(W))
        W <- matrix(rnorm(L, sd=.5), L, NCOL(y))

    PF("NWT.LMM.DV1.BEGIN = \n")
    print(round(lmm.dv1(W, K, y)[[1]], 4))

    time0 <- proc.time()
    for(i in seq(100))
    {
        g <- lmm.dv1(W, K, y)[[1]]
        if(max(abs(g)) < 1e-6)
            break
        H <- lmm.dv2(W, K, y)[[1]]

        ## update: u = -g H^{-1} => H u = -g
        u <- solve(H, -g)
        if(max(abs(u)) < 1e-6)
            break
        W <- W + u
    }
    time1 <- proc.time()

    ## report
    PF("NWT.LMM.DV2.END =\n")
    print(round(H, 4))

    dv1 <- lmm.dv1(W, K, y)
    PF("NWT.LMM.DV1.END =\n")
    print(round(g, 4))

    ## make predictions
    prd <- knl.prd(y, K, W)
    
    ## timing
    rtm <- DF(key='rtm', val=unname((time1 - time0)['elapsed']))

    list(rpt=rbind(rtm, prd), par=W)
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
