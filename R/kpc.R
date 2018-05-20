## Kernel Perceptron (KPC)

#' @title Gradient of LSE loss wrt. Deep Kernel Net parameters.
#' @param luy vector J weights that mix J inner kernels of {U}, the latent data, to a single
#' covariance of {Y}, the output;
#' @param lxu matrix L by M weights that mix L basic kernels of {X}, the (unknown) input, to
#' M covariances of {U[, 1..M]}, the M lantent features;
#' @param phi scalar size of the universal noisy;
#' @param y vector the observed, noisy responses;
#' @param bkn list of basic kernel matrices;
#' @param ikn list of inner kernel names or constructors;
#' @param S integer the number of N-MVN samples to draw to estimate mean wrt. the latent U.
up0 <- function(par, ctx, S=1e3, dvt=TRUE, wfn=softplus, wdv=logistic, ...)
{
    ## training data and network structure.
    bkn <- ctx$bkn
    y <- ctx$y

    ## dimensions
    L <- length(bkn)                    # number of inner kernels
    N <- NROW(y)                        # number of individuals
    
    ## parameters
    phi <- par[L+1] %>% unname
    txy <- par[1:L] %>% unname
    PHI <- wfn(phi)
    TXY <- wfn(txy)
    
    ## noise covariance
    IMX <- diag(1.0, N)

    ## L basic kernels organized in an N-N-L array
    BKN <- vapply(c(bkn, list(IMX)), unname, IMX)
    
    ## Y's covariance and accuracy (inverse covariance)
    YCV <- BKN %>% .dim(N^2, L+1) %>% .mbm(c(TXY, PHI)) %>% .dim(N, N)
    YCL <- chol(YCV)
    YAC <- chol2inv(YCL)
    YAY <- drop(YAC %*% y)
    
    ## prediction & loss
    yht <- (YCV - diag(PHI, N)) %*% YAY

    ## objective
    luy <- .5 * sum(crossprod(y, YAC) * y) + sum(log(diag(YCL))) + .5 * N * log(2*pi)
    
    ctx <- within(ctx, {yht=yht; luy=luy})
    if(!dvt)
        return(ctx)

    ## derivatives of L(Y|U) wrt. tau_{L}
    vapply(1:L, function(l)
    {
        -.5 * sum((tcrossprod(YAY) - YAC) * BKN[, , l])
    },
    numeric(1)) * wdv(txy) -> dvt.lxy
    
    ## derivatives of L(Y|U) wrt. phi
    dvt.phi <- -.5 * sum(YAY^2, - diag(YAC)) * wdv(phi)
    
    ## derivatives
    dvt <- list(lxy=dvt.lxy, phi=dvt.phi) %>% unlist %>% {./N}
    
    ## update context and return
    ctx <- within(ctx, {dvt <- dvt; par <- par})
    ctx
}

#' @title predict noise free outcome
#' most of the procedure are the same as update, but since no gradient calculation is
#' required, this function runs faster.
pd0 <- function(par, ctx, wfn=sq, ...)
{
    up0(par, ctx, wfn, dvt=FALSE, ...)$yht
}

#' draw mini-batches from the full kernel matrix
.bat <- function(ctx, bsz=47, itr=0)
{
    N <- length(ctx$y)
    ## not by batch
    if(is.null(bsz) || bsz >= N)
        bsz <- N

    ## advance epoch index
    ep <- (itr * bsz^2) / N^2           # epoch count
    bt <- (itr * bsz^2) %% N^2 / bsz^2  # batch count
    ## ep <- (itr * bsz) / N
    ## bt <- (itr * bsz) %% N /bsz
    ## randomly pick out {bsz} individuals
    if(bsz < N)
    {
        idx <- sample.int(N, bsz)
        ctx <- within(ctx,
        {
            bkn <- lapply(bkn, `[`, idx, idx)
            y <- y[idx]
        })
    }
    ctx <- within(ctx, {ep <- ep; bt <- bt})
    ctx
}

    
#' @title gradient descent for kernel deep net
#' @param par vector initial values of parameters;
#' @param ctx list context of the training task, such as input features and lables;
#' @param usz integer number of sampling to calculate the losses and derivatives;
#' @param bsz integer batch size of stochastic gradient descent;
#' @param mmt numeric momment of stochastic gradient descent, between 0 and 1;
#' @param lr0 numeric the initial learning rate;
#' @param wtm numeric (in ...) walltime in hours
#' @param esz integer (in ...) number of samples to draw in for error estimation;
#' also collect training and evaluation statistics such as L1 and L2 error.
gd0 <- function(par, ctx, usz=1e3, bsz=30, lr0=1e-3, ...)
{
    ## list of history
    hst <- list()

    ## dots
    dot <- list(...)
    max.itr <- dot$max.itr %||% 100
    wtm <- dot$wtm %||% 5e-2
    tol <- dot$tol %||% 1e-6
    mmt <- dot$mmt %||% 0
    esz <- dot$esz %||% usz

    ## evaluation dataset
    y <- ctx$y
    evl <- ctx$evl

    ## initial parameters
    lrt <- lr0
    dvt <- 0 * par
    dff <- 0
    i <- 1
    time.start <- Sys.time()
    while(TRUE)
    {
        stt <- .bat(ctx, bsz, i-1)
        ep <- stt$ep
        bt <- stt$bt
        ## ep <- i - 1
        ## bt <- 0
        
        tm0 <- Sys.time()
        stt <- up0(par, stt, usz, ...)
        sec <- (Sys.time() - tm0) %>% {units(.) <- 'secs';.} %>% as.numeric
        dvt <- (mmt * dvt) + (1 - mmt) * stt$dvt
        
        ## update statistics upon new epoch
        if(bt < 1)
        {
            ## for developing data
            if(!is.null(bsz) && bsz < NROW(y))
                yht <- pd0(stt$par, ctx, esz=esz, ...)
            else
                yht <- stt$yht
            dve <- sum((ctx$y - yht)^2) / length(yht)
            dvc <- cor(ctx$y, yht)
            luy <- stt$luy
        }

        ## record each iteration
        hst[[i]] <-  list(
            i=i, ep=ep, bt=bt, sec=sec,
            par=par, dvt=dvt, lrt=lrt,
            dve=dve, dvc=dvc)

        ## update learning rate
        if(i == 1)
        {
            lrt <- getLearningRate(lr = lrt, type = "Specified")
            ## header of tracks
            hdr <- c(
                eph=sprintf("%4s", 'eph'),     # epoch
                bat=sprintf("%4s", 'bat'),     # batch
                luy=sprintf("%7s", 'luy'),
                dve=sprintf('%7s', 'dve'),     # developing erro
                dvc=sprintf('%4s', 'dvc'),     # developing corr
                dvt=sprintf('%9s', 'xtm.dvt'), # extreme derivative
                mmt=sprintf('%3s', 'mmt'),     # mommentum
                lrt=sprintf('%9s', 'lrt'),     # learning rate
                phi=sprintf('%9s', 'par.phi'), # phi
                sec=sprintf('%4s', 'sec'))     # seconds used
            cat(paste(hdr, collapse=' '), '\n', sep='')
        }

        ## message tracks
        xmt.dvt <- dvt[which.max(abs(dvt))]
        msg <- c(
            eph=sprintf("%04X", as.integer(ep)), # epoch
            bat=sprintf("%04X", as.integer(bt)), # batch
            luy=sprintf('%7.2f', luy),           # error
            dve=sprintf('%7.2f', dve),           # error
            dvc=sprintf('%4.2f', dvc),           # corr(y, yhat)
            dvt=sprintf('%9.1e', xmt.dvt),       # extreme derivative
            mmt=sprintf('%3.1f', mmt),           # mommentum
            lrt=sprintf('%9.1e', lrt),           # learning rate
            phi=sprintf('%9.3f', par['phi']),    # phi
            sec=sprintf('%4.1f', sec),           # seconds used
            dff=sprintf('%9.1e', dff))           # 
        cat(paste(msg, collapse=' '), '\n', sep='')

        if(i > max.itr)
        {
            cat('SGD: reaching max iter:', max.itr, '\n')
            break
        }
        if((Sys.time() - time.start) %T>% {units(.) <- 'hours'} %>% as.numeric > wtm)
        {
            cat('SGD: reaching walltime:', wtm, 'hour(s)\n')
            break
        }

        ## greadient descent here:
        par <- par - lrt * dvt
        dff <- sum(abs(unlist(par) - unlist(hst[[i]]['par'])))
        if(dff < tol)
            break
        i <- i + 1
    }

    ## rearrange history
    hst <- hst %>% lapply(unlist) %>% {do.call(rbind, .)} %>% as.data.frame %>% as.tbl
    hst <- hst %>% mutate(sec=cumsum(sec) - first(sec))
    ## drop paramters, gradient, and fractional epoches
    hst <- hst %>% filter(bt == 0) %>% select(-bt)
    hst <- hst %>% select(-starts_with('par.'), -starts_with('dvt.'))
    
    ## return the history and new parameters
    ret <- within(stt,
    {
        dve <- dve
        dvc <- dvc
        dvt <- dvt
        par <- par
        hst <- hst
    })
    ret
}
