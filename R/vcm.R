## variance component model (VCM)
vcm.dv1 <- function(W, K, Y, ...)
{
    M <- NCOL(Y)                        # dim of Y, upper units
    L <- NROW(W)                        # dim of K, lower kernels
    W <- as.matrix(exp(W))
    C <- cmb(K, W)                      # cov: L x M cov matrices
    Y <- as.matrix(Y)

    ## 1st derivative
    dv1 <- sapply(seq.int(M), function(m)
    {
        ## a_i
        A <- chol2inv(chol(C[[m]]))

        ## alpha_i = a_i %*% y_i
        YAY <- A %*% Y[, m]

        ## tmp_i = alpha_i alpha_i' - a_i
        TMP <- tcrossprod(YAY) - A
        sapply(seq.int(L), function(l) -.5 * sum(TMP * K[[l]]) * W[l, m])
    })
    dv1
}

cpp.dv1 <- function(W, K, Y, ...)
{
    W <- as.matrix(W)
    Y <- as.matrix(Y)
    if(!is.list(K))
        K <- list(K)
    .Call('vcm_dv1', W, K, Y)
}


vcm.dv2 <- function(W, K, Y, ...)
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
        TMP2 <- 2 * YAY2[[m]] - A[[m]]
        TMP3 <- YAY2[[m]] - A[[m]]
        for(r in seq.int(1, L))
        {
            for(s in seq.int(r, L))
            {
                dmx[r, s] <- .5 * W[r, m] * sum(A[[m]] %*% K[[r]] * TMP2 %*% K[[s]]) * W[s, m]
                dmx[s, r] <- dmx[r, s]
            }
            dmx[r, r] <- dmx[r, r] -.5 * sum(TMP3 * K[[r]]) * W[r, m]
        }
        dmx
    })
    dv2
}

vcm.fsi <- function(W, K, Y, ...)
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

#' Variance component model prediction
#' 
#' @param y a vector of response variable
#' @param K a list of covariance kernels
#' @param W a vector of variance components
#' @param ln the VCs are logged? if true, exp(W) is used instead of W
#' @param rt return a table? if false, a named vector is returned
vpd <- function(y, K, W, rt=1, ...)
{
    y <- unname(y)
    ## is W logged
    N <- NROW(y)

    ## prepand noisy kernel
    C <- c(list(eps=diag(N)), cv=K)
    
    ## make predictions
    v <- unname(cmb(C, W)[[1]])
    ## a <- solve(v)
    alpha <- solve(v, y)                # V^{-1}y
    f <- v - diag(W[1], length(y))      # W[1] is PHI

    ## prediction 1: conditional Gaussian
    h <- f %*% alpha

    mse <- mean((y - h)^2)
    ## gurad against zero-standard deviation
    cyh <- tryCatch(cor(y, h), warning=function(w) 0, error=function(e) NA)

    ## negative likelihood
    ld <- with(.Internal(det_ge_real(v, TRUE)), sign * modulus)
    nlk <- .5 / N * (sum(alpha * y) + ld + N * log(2*pi))
    
    ## prediction 2: leave one out CV
    ## h <- y - a %*% y / diag(a)
    ## loo <- mean((y - h)^2)

    ## return
    rpt <- c(mse=mse, nlk=nlk, cyh=cyh, ssz=N)
    if(rt == 1)
        rpt <- DF(key=names(rpt), val=rpt, row.names=names(rpt))
    if(rt == 2)
        rpt <- DF(t(rpt))
    rpt
}

rop.vcm <- function(y, K, W=NULL, cpp=TRUE, ...)
{
    N <- NROW(y)
    Q <- NCOL(y)
    C <- c(list(eps=diag(N)), cv=K)
    L <- length(C)
    if(is.null(W))
        W <- matrix(rnorm(L * Q), L, Q)

    obj <- function(x)
    {
        v <- cmb(C, exp(x))[[1]]
        u <- try(chol(v))
        if(inherits(u, 'try-error'))
        {
            alpha <- solve(v, y)                # V^{-1}y
            v.det <- with(.Internal(det_ge_real(v, TRUE)), sign * modulus)
        }
        else
        {
            alpha <- backsolve(u, forwardsolve(u, y, upper.tri=TRUE, transpose=TRUE))
            v.det <- 2 * sum(log(diag(u)))
        }
        .5 * (sum(alpha * y) + v.det + N * log(2*pi))
    }
    if(cpp)
        grd <- function(x) cpp.dv1(x, C, y)[, 1]
    else
        grd <- function(x) vcm.dv1(x, C, y)[, 1]
    hsn <- function(x) vcm.dv2(x, C, y)[[1]]
    fsi <- function(x) vcm.fsi(x, C, y)[[1]]

    ## using R's optimizer
    ## library(numDeriv)
    ## print(list(grd.num=grad(obj, W), grd.fun=grd(W)))
    ## print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))
    time0 <- proc.time()
    ret <- optim(W, obj, grd, method="L-BFGS-B", control=list(trace=0))
    time1 <- proc.time()
    W <- ret$par
    ## print(list(grd.num=grad(obj, W), grd.fun=grd(W)))
    ## print(list(hsn.num=hessian(obj, W), hsn.fun=hsn(W), fsn.fun=fsi(W)))

    ## make predictions
    prd <- vpd(y, K, exp(W))

    ## timing
    rtm <- DF(key='rtm', val=unname((time1 - time0)['elapsed']))

    ret$rpt <- rbind(rtm=rtm, prd)
    ret$par <- exp(W)
    ret
}

nwt.vcm <- function(y, K, W=NULL, ...)
{
    . <- list(...)
    wep <- .$wep %||% 100
    mmt <- .$mmt %||% 1.0
    vrb <- .$vrb %||% 1

    N <- NROW(y)
    C <- c(list(eps=diag(N)), cv=K)
    L <- length(C)
    Q <- NCOL(y)
    if(is.null(W))
        W <- matrix(rnorm(L * Q), L, Q)

    PF("NWT.VCM.DV1.BEGIN = \n")
    time0 <- proc.time()
    for(i in seq.int(wep))
    {
        g <- cpp.dv1(W, C, y)[, 1, drop=FALSE]
        if(max(abs(g)) < 1e-6)
            break

        H <- vcm.fsi(W, C, y)[[1]]
        ## print(list(itr=i, dv2=H, dv1=g, par=exp(W)))

        ## update: u = -g H^{-1} => H u = -g
        u <- try(solve(H, -g))
        if(inherits(u, 'try-error'))
            break
        if(max(abs(u)) < 1e-6)
            break
        W <- W + u
    }
    time1 <- proc.time()
    PF("NWT.VCM.DV2.END =\n")
    
    ## make predictions
    prd <- vpd(y, K, exp(W))
    
    ## timing
    rtm <- DF(key='rtm', val=unname((time1 - time0)['elapsed']))

    list(rpt=rbind(rtm, prd), par=exp(W))
}

vcm <- function(y, v, e)
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
    h <- y - a %*% y / diag(a)
    loo <- mean((y - h)^2)

    ## negative log likelihood
    nlk <- nlk(y, v, u, a)

    DF(key=c('mse', 'loo', 'nlk', 'cyh'), val=c(mse, loo, nlk, cyh))
}

#' the null model prediction
nul <- function(y)
{
    N <- NROW(y)

    ## the distribution of y
    m <- mean(y)
    v <- diag(var(y), N)

    ## mean square error
    mse <- mean((y - m)^2)

    ## negative log likelihood
    nlk <- nlk(y, v) / N

    ## leave one out
    h <- (sum(y) - y) / (N - 1)
    loo <- mean((y - h)^2)

    ## report
    ret <- DF(key=c('mse', 'nlk', 'loo', 'cyh'), val=c(mse, nlk, loo, cyh))
    rownames(ret) <- ret$key
    ret
}

nul.vcm <- function(y, eps=NULL, ...)
{
    N <- NROW(y)
    mse <- sum(y^2) / N

    if(is.null(eps))
        eps <- mse
    v <- diag(eps, N)
    alpha <- solve(v, y)                # v^{-1}y or y/e = y / sum(y^2/N) = N * y/ sum(y^2)

    ldt <- with(.Internal(det_ge_real(v, TRUE)), sign * modulus)
    nlk <- .5 / N * (sum(alpha * y) + ldt + N * log(2*pi))
    
    rpt <- DF(key=c('mse', 'nlk'), val=c(mse, nlk))
    rownames(rpt) <- rpt$key

    ## return
    list(rpt=rpt, par=c(eps=eps))
}