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
#' @param v list of covariance kernels
#' @param w vector of parameters (beta, and sigma^2)
#' @param x matrix of covariate
#' @param ... additional arguments
vpd <- function(y, v=NULL, x=NULL, w, ...)
{
    dot <- list(...)
    y <- unname(y)
    N <- NROW(y)
    w[is.na(w)] <- 0
    
    ## kernels, prepended with noise
    v <- c(list(EPS=diag(N)), v)
    v <- v[intersect(names(v), names(w))]
    K <- length(v)

    ## fix effects, prepended with intercept
    x <- cbind(X00=rep(1, N), x)
    M <- ncol(x)

    ## variance components, and fixed coeficients
    w <- unlist(w)
    vcs <- w[intersect(names(w), names(v))]
    fix <- w[intersect(names(w), colnames(x))]
    
    ## fixed effect, and residual
    xb <- if(M > 0) x %*% fix else 0    # x beta
    e <- y - xb                         # fix effect residual
    SST <- sum(e^2) / (N - M)           # sum square total

    ## heritability
    hsq <- unname(1 - vcs[1] / SST)      # hsq 2
    
    ## predicted random effects
    V <- unname(cmb(v, vcs, TRUE))
    A <- try(chol2inv(chol(V)), silent=TRUE)
    if(inherits(A, 'try-error'))
    {
        A <- MASS::ginv(V)
    }
    Ae <- A %*% e
    
    zub <- (V - diag(vcs[1], N)) %*% Ae   # use BLUP
    zul <- e - Ae / diag(A)               # LOO
    yhb <- zub + xb
    yhl <- zul + xb

    zeb <- mean((e - zub)^2)     # e MSE 1, BLUP
    zel <- mean((e - zul)^2)     # e MSE 2, LOOV
    yeb <- mean((y - yhb)^2)     # y MSE 1, BLUP
    yel <- mean((y - yhl)^2)     # y MSE 2, LOOV

    ## correlation between truth and prediction
    ycb <- if(abs(sd(yhb)) < 1e-8) 0 else cor(y, yhb)
    ycl <- if(abs(sd(yhl)) < 1e-8) 0 else cor(y, yhl)
    zcb <- if(abs(sd(zub)) < 1e-8) 0 else cor(e, zub)
    zcl <- if(abs(sd(zul)) < 1e-8) 0 else cor(e, zul)

    ## negative likelihood
    ldt <- with(determinant(V), modulus * sign) / N
    eae <- sum(Ae * e) / N    # e^T V^{-1} e
    ## nlk <- .5 * (eae + ldt + log(2 * pi))
    nlk <- eae + ldt
    attributes(nlk) <- NULL
    
    ## return
    rpt <- data.frame(N=N, hsq=hsq, SST=SST, nlk=nlk,
                      zeb=zeb, zel=zel, yeb=yeb, yel=yel,
                      ycb=ycb, ycl=ycl, zcb=zcb, zcl=zcl)

    ## 1st half predict 2nd half
    i <- sample(c(TRUE, FALSE), N, replace=TRUE)
    j <- !i
    f <- e; f[i] <- NA; zui <- kpd(f, V)
    f <- e; f[j] <- NA; zuj <- kpd(f, V)
    zuh <- numeric(N)
    zuh[i] <- zui
    zuh[j] <- zuj
    yhh <- zuh + xb
    rpt$yeh <- mean((y - yhh)^2)
    rpt$ych <- if(abs(sd(yhh)) < 1e-8) 0 else drop(cor(y, yhh))
    rpt
}

#' Known outcome predict unknown
kpd <- function(e, V)
{
    j <- is.na(e)                       # known
    i <- !j                             # unknown
    i <- which(i)
    j <- which(j)
    
    z <- drop(V[j, i] %*% solve(V[i, i], e[i]))
    z
}

loo <- function(y, v=NULL, x=NULL, w, ...)
{
    dot <- list(...)
    y <- unname(y)
    N <- NROW(y)
    Q <- length(w)
    
    ## kernels, prepended with noise
    v <- c(list(EPS=diag(N)), v)
    v <- v[intersect(names(v), names(w))]
    K <- length(v)

    ## fix effects, prepended with intercept
    x <- cbind(X00=rep(1, N), x)
    M <- ncol(x)

    ## variance components, and fixed coeficients
    vcs <- w[intersect(names(w), names(v))]
    fix <- w[intersect(names(w), colnames(x))]
    
    ## fixed effect, and residual
    xb <- if(M > 0) x %*% fix else 0    # x beta
    e <- y - xb                         # fix effect residual
    SST <- sum(e^2) / (N - M)           # sum square total

    ## heritability
    hsq <- unname(1 - vcs[1] / SST)      # hsq 2
    ## h2c <- unname(1 - vcs[1] / sum(vcs)) # hsq 3
    
    ## predicted random effects
    V <- unname(cmb(v, vcs, TRUE))
    A <- try(chol2inv(chol(V)), silent=TRUE)
    if(inherits(A, 'try-error'))
       A <- ginv(V)
    Ae <- A %*% e
    
    zul <- e - Ae / diag(A)               # LOO
    yhl <- zul + xb
    yhl
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

    rpt <- vpd(y, w=par, rt=0)
    list(rpt=rpt, par=par, rtm=NA)
}
