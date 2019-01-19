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
#' @param y vector of response variable
#' @param W vector of parameters (beta, and sigma^2)
#' @param K list of covariance kernels
#' @param X matrix of covariate
#' @param rt return a table? if false, a named vector is returned
vpd <- function(W, y, K=NULL, X=NULL, rt=1, ...)
{
    dot <- list(...)
    y <- unname(y)
    N <- NROW(y)
    Q <- length(W)

    ## kernels, prepended with noise
    K <- c(list(EPS=diag(N)), K)
    L <- length(K)
    M <- if(is.null(X)) 0 else NCOL(X)

    ## variance components, and fix effect coeficients
    vcs <- W[seq(Q - L + 1, Q)]
    fix <- W[seq(1, Q - L)]
    
    ## matrix for fixed effect, prepended with intercept
    if(length(fix) == M + 1)
    {
        X <- cbind(X00=rep(1, N), X)
        M <- M + 1
    }

    xb <- if(M > 0) X %*% fix else 0    # x beta
    e <- y - xb                         # fix effect residual
    SST <- sum(e^2) / (N - 1)           # sum square total

    ## heritability
    h2b <- unname(1 - vcs[1] / SST)      # hsq 2
    h2c <- unname(1 - vcs[1] / sum(vcs)) # hsq 3
    
    ## predicted random effects
    V <- unname(cmb(K, vcs, drop=TRUE))
    A <- solve(V)                       # V^{-1}
    Ae <- A %*% e
    
    zu1 <- (V - diag(vcs[1], N)) %*% Ae   # use BLUP
    zu2 <- e - Ae / diag(A)               # LOO
    yh1 <- zu1 + xb
    yh2 <- zu2 + xb

    er1 <- mean((e - zu1)^2)            # e MSE 1, BLUP
    er2 <- mean((e - zu2)^2)            # e MSE 2, LOOV
    ey1 <- mean((y - yh1)^2)            # y MSE 1, BLUP
    ey2 <- mean((y - yh2)^2)            # y MSE 2, LOOV
    
    ## correlation between truth and prediction
    cy1 <- cor(y, yh1)
    cy2 <- cor(y, yh2)
    cr1 <- if(abs(sd(zu1)) < 1e-8) 0 else cor(e, zu1)
    cr2 <- if(abs(sd(zu2)) < 1e-8) 0 else cor(e, zu2)

    ## slops between truth and prediction
    sr1 <- if(abs(sd(zu1)) < 1e-8) 0 else unname(coef(lm(e ~ zu1))[2])
    sr2 <- if(abs(sd(zu2)) < 1e-8) 0 else unname(coef(lm(e ~ zu2))[2])
    sy1 <- unname(coef(lm(y ~ yh1))[2])
    sy2 <- unname(coef(lm(y ~ yh2))[2])
    
    ## negative likelihood
    ldt <- with(.Internal(det_ge_real(V, TRUE)), sign * modulus) / N
    eae <- sum(Ae * e) / N    # e^T V^{-1} e
    ## nlk <- .5 * (eae + ldt + log(2 * pi))
    nlk <- eae + ldt
    
    ## return
    rpt <- c(h2b=h2b, h2c=h2c, SST=SST, nlk=nlk,
             er1=er1, er2=er2, ey1=ey1, ey2=ey2,
             cy1=cy1, cy2=cy2, cr1=cr1, cr2=cr2,
             sr1=sr1, sr2=sr2, sy1=sy1, sy2=sy2, N=N)

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
    C <- c(list(EPS=diag(N)), cv=K)
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
    prd <- vpd(exp(W), y, K)

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
    C <- c(list(EPS=diag(N)), cv=K)
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
    prd <- vpd(exp(W), y, K)
    
    ## timing
    rtm <- DF(key='rtm', val=unname((time1 - time0)['elapsed']))

    list(rpt=rbind(rtm, prd), par=exp(W))
}

#' Null Variance Component Model
#'
#' @param y response
#' @param X design matrix of fix effect, without intercept.
#'
#' @return null model fit.
nul.vcm <- function(y, X=NULL, ...)
{
    N <- NROW(y)
    m0 <- if(is.null(X)) lm(y ~ 1) else lm(y ~ X)

    ## parameters
    par <- coef(m0)
    names(par) <- if(is.null(X)) 'X00' else c('X00', names(X))
    par <- c(par, EPS=sum(residuals(m0)^2) / (N - 1))

    rpt <- vpd(par, y, rt=0)
    list(rpt=rpt, par=par, rtm=NA)
}
