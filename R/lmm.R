## LMM or pure variance component
lmm.dv1 <- function(W, K, Y, ...)
{
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- matrix(exp(W))

    C <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(C, chol)                # decompotion
    A <- lapply(H, chol2inv)            # acc: a_i, (i = 1 .. M)

    Y <- as.matrix(Y)

    ## alpha_i = a_i %*% y_i
    YAY1 <- lapply(1:M, function(m) A[[m]] %*% Y[, m])

    ## tmp_i = alpha_i alpha_i' - a_i
    TMP <- lapply(1:M, function(m) tcrossprod(YAY1[[m]]) - A[[m]])

    ## 1st derivative
    dv1 <- lapply(seq.int(M), function(m)
    {
        sapply(seq.int(L), function(l)
        {
            -.5 * sum(TMP[[m]] * K[[l]]) * W[l, m]
        })
    })
    dv1
}

lmm.dv2 <- function(W, K, Y, ...)
{
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- matrix(exp(W))

    C <- cmb(K, W)                      # cov: L x M cov matrices
    H <- lapply(C, chol)                # decompotion
    A <- lapply(H, chol2inv)            # acc: a_i, (i = 1 .. M)

    Y <- as.matrix(Y)

    ## ## alpha_i = a_i y_i
    YAY1 <- lapply(1:M, function(m) A[[m]] %*% Y[, m])

    ## beta_i = (a_i y_i) (a_i y_i)' = alpha_i alpha_i'
    YAY2 <- lapply(YAY1, tcrossprod)
    
    ## container of second derivative
    dv2 <- lapply(1:M, function(m)
    {
        dmx <- matrix(.0, L, L)
        for(r in seq.int(1, L))
        {
            for(s in seq.int(r, L))
            {
                TMP2 <- 2 * YAY2[[m]] - A[[m]]
                dmx[r, s] <- .5 * W[r, m] * sum(A[[m]] %*% K[[r]] * TMP2 %*% K[[s]]) * W[s, m]
                dmx[s, r] <- dmx[r, s]
            }
            TMP2 <- YAY2[[m]] - A[[m]]
            dmx[r, r] <- dmx[r, r] -.5 * sum(TMP2 * K[[r]]) * W[r, m]
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
    W <- matrix(exp(W))

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

#' @param y a vector of response variable
#' @param K a list of covariance kernels
#' @param W a vector of variance components
#' @param ln the VCs are logged? if true, exp(W) is used instead of W
#' @param rt return a table? if false, a named vector is returned
knl.prd <- function(y, K, W, ln=1, rt=1, ...)
{
    ## is W logged
    if(ln)
        W <- as.matrix(exp(W))
    N <- NROW(y)

    ## prepand noisy kernel
    C <- c(list(eps=diag(N)), cv=K)
    
    ## make predictions
    v <- cmb(C, W)[[1]]
    f <- v - diag(W[1], length(y))      # W[1] is PHI

    ## prediction 1: conditional Gaussian
    y <- unname(drop(as.vector(y)))
    h <- unname(drop(f %*% (solve(v) %*% y)))

    mse <- mean((y - h)^2)
    ## gurad against zero-standard deviation
    cyh <- tryCatch(cor(y, h), warning=function(w) 0, error=function(e) NA)
    nlk <- nlk(y, v)
    ## return
    rpt <- c(mse=mse, nlk=nlk, cyh=cyh, ssz=N)
    if(rt == 1)
        rpt <- DF(key=names(rpt), val=rpt)
    if(rt == 2)
        rpt <- DF(t(rpt))
    rpt
}

rop.lmm <- function(y, K, W=NULL, ...)
{
    N <- NROW(y)
    Q <- NCOL(y)
    C <- c(list(eps=diag(N)), cv=K)
    L <- length(C)
    if(is.null(W))
        W <- matrix(rchisq(L * Q, df=1) * .1, L, Q)
    W <- log(W)

    obj <- function(x) nlk(y, cmb(C, exp(x))[[1]])
    grd <- function(x) lmm.dv1(x, C, y)[[1]]
    hsn <- function(x) lmm.dv2(x, C, y)[[1]]
    fsi <- function(x) lmm.fsi(x, C, y)[[1]]

    ## using R's optimizer
    ## library(numDeriv)
    ## print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))
    ## dv1 <- lmm.dv1(W, C, y)[[1]]        # gradient before
    ## PF("DV1.MAX = %9.3f\n", dv1[which.max(abs(dv1))])
    time0 <- proc.time()
    ret <- optim(W, obj, grd, method="L-BFGS-B", control=list(trace=0))
    time1 <- proc.time()
    W <- ret$par
    ## print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))
    ## dv1 <- lmm.dv1(W, C, y)[[1]]        # gradient after
    ## PF("DV1.MAX = %9.3f\n", dv1[which.max(abs(dv1))])

    ## make predictions
    prd <- knl.prd(y, K, exp(W), ln=0)

    ## timing
    rtm <- DF(key='rtm', val=unname((time1 - time0)['elapsed']))

    ret$rpt <- rbind(rtm, prd)
    ret$par <- exp(W)
    ret
}

nwt.lmm <- function(y, K, W=NULL, ...)
{
    . <- list(...)
    wep <- .$wep %||% 100
    mmt <- .$mmt %||% 1.0
    vrb <- .$vrb %||% 1

    N <- NROW(y)
    C <- c(list(eps=diag(N)), cv=K)
    L <- length(C)
    if(is.null(W))
        W <- matrix(rnorm(L, sd=.5), L, NCOL(y))

    PF("NWT.LMM.DV1.BEGIN = \n")
    print(round(lmm.dv1(W, C, y)[[1]], 4))

    time0 <- proc.time()
    for(i in seq.int(wep))
    {
        g <- lmm.dv1(W, C, y)[[1]]
        if(max(abs(g)) < 1e-5)
            break

        H <- lmm.fsi(W, C, y)[[1]]
        ## print(list(itr=i, hsn=H))

        ## update: u = -g H^{-1} => H u = -g
        u <- try(solve(H, -g), TRUE)
        if(inherits(u, 'try-error'))
            break
        if(max(abs(u)) < 1e-5)
            break
        W <- W + u
    }
    time1 <- proc.time()

    ## report
    print(list(dv2=round(H, 3), dv1=round(g, 3), par=W))
    PF("NWT.LMM.DV2.END =\n")
    
    ## make predictions
    prd <- knl.prd(y, K, W)
    
    ## timing
    rtm <- DF(key='rtm', val=unname((time1 - time0)['elapsed']))

    list(rpt=rbind(rtm, prd), par=W)
}
