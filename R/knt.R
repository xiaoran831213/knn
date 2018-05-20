## linear combination of kernels
cmb <- function(k, w)
{
    if(is.matrix(k))
        k <- list(k)
    w <- as.matrix(w)
    M <- ncol(w)
    L <- nrow(w)
    lapply(1:M, function(m) Reduce(`+`, mapply(`*`, k, w[, m], SIMPLIFY=FALSE)))
}

mky <- function(ctx, rep=1000, ...)
{
    for(. in names(ctx)) assign(., ctx[[.]])
    N <- NROW(ukn[[1]])                        # number of individuals
    TUY <- exp(par$tuy)
    TXU <- exp(par$txu)
    ## L basic kernels organized in an N-N-L array
    ucv <- cmb(ukn, TXU)
    ucl <- lapply(ucv, chol)
    uac <- lapply(ucl, chol2inv)

    ## calculate mean ycv over umx
    ycv <- matrix(.0, N, N)
    for(.i in seq.int(rep))
    {
        U <- vapply(1:M, function(m) mvn(1, ucl[[m]], 1), rep(.0, N))
        ikn <- findKernel(U, ykn)
        ycv <- ycv + cmb(ikn, TUY)[[1]]
    }
    ycv <- ycv / rep
    ycv
}

#' @title Gradient of LSE loss wrt. Deep Kernel Net parameters.
#' @param luy vector J weights that mix J inner kernels of {U}, the latent data, to a single
#' covariance of {Y}, the output;
#' @param lxu matrix L by M weights that mix L basic kernels of {X}, the (unknown) input, to
#' M covariances of {U[, 1..M]}, the M lantent features;
#' @param phi scalar size of the universal noisy;
#' @param y vector the observed, noisy responses;
#' @param ukn list of basic kernel matrices;
#' @param ikn list of inner kernel names or constructors;
#' @param S integer the number of N-MVN samples to draw to estimate mean wrt. the latent U.
up3 <- function(.cx, .dv=2, ...)
{
    for(. in names(.cx)) assign(., .cx[[.]])
    L <- length(ukn)                    # number of basic kernels
    J <- length(ykn)                    # number of inner kernels
    N <- NROW(y)                        # number of individuals
    TXU <- exp(par$txu)
    TUY <- exp(par$tuy)
    ## L basic kernels organized in an N-N-L array
    ucv <- cmb(ukn, TXU)
    ucl <- lapply(ucv, chol)
    uac <- lapply(ucl, chol2inv)
    
    if(!exists('U', inherit=FALSE))     # the very first hidden data
    {
        U <- vapply(1:M, function(m) mvn(1, ucl[[m]], 1), y)
        U <- rmf.matrix(U)
    }
    ikn <- findKernel(U, ykn)
    ycv <- cmb(ikn, TUY)[[1]]
    ycl <- chol(ycv)
    yac <- chol2inv(ycl)
    yay <- drop(yac %*% y)
    
    dvt <- lapply(par, `*`, 0)
    if(.dv == 1 || .dv > 3)             # advance upper level
    {
        dvt$tuy[,] <- vapply(1:J, function(j)
        {
            -.5 * sum((tcrossprod(yay) - yac) * ikn[[j]])
        },
        .0) * exp(par$tuy)
    }
    if(.dv == 2)                        # propose hidden data
    {
        for(.i in 1:10)
        {
            ## Propose a walkaway from U
            ## .mx <- vapply(1:M, function(m) mvn(1L, ucl[[m]], 1L), y) * .02 + U
            .rt <- rrot(N, 0, .5)
            .mx <- .rt %*% U
            .cv <- cmb(findKernel(.mx, ykn), TUY)[[1]]
            .cl <- chol(.cv)
            .ac <- chol2inv(.cl)
            
            ## log[f(u1, y) / f(u, y)], or log[f(u1|y)/f(u|y)]
            .pu <- sapply(1:M, function(m)
            {
                -.5 * xvy(.mx[, m] - U[, m], uac[[m]], .mx[, m] + U[, m])
            }) %>% sum
            .py <- -.5 * xvy(y, .ac - yac) - sum(log(diag(.cl)) - log(diag(ycl)))
            .pp <- .pu + .py            #

            ## if(luy < lxu && .pu > max(0, .py) || luy > lxu && .py > max(0, .pu))
            if(.pp > 0)
            {
                U <- .mx
                ycv <- .cv
                ycl <- .cl
            }
        }
    }
    uau <- vapply(1:M, function(m) uac[[m]] %*% U[, m], y)
    if(.dv == 3 || .dv > 3)             # advance lower level
    {
        dvt.lxu <- matrix(.0, L, M)
        for(.m in seq.int(M))
        {
            .tmp <- tcrossprod(uau[, .m]) - uac[[.m]]
            for(.l in seq.int(L))
            {
                dvt.lxu[.l, .m] <- -.5 * sum(.tmp * ukn[[.l]])
            }
        }
        dvt$txu <- dvt.lxu * exp(par$txu)
    }
    dvt <- lapply(dvt, `/`, N)
    
    ## predictions (hats)
    .uh <- sapply(1:M, function(m)
    {
        (ucv[[m]] - diag(TXU[1, m], N)) %*% (uac[[m]] %*% U[, m])
    })
    .uh <- .uh / sqrt(colSums(.uh^2))
    .vh <- cmb(findKernel(.uh, ykn), TUY)[[1]]
    .cl <- chol(.vh)
    .ac <- chol2inv(.cl)
    ## .ay <- drop(.ac %*% y)
    ## yht <- (.vh - diag(TUY[1], N)) %*% .ay
    ## ylo <- y - .ay / diag(.ac)
    ## training statistics
    lxu <- -sum(sapply(1:M, function(m) .mlk(U[, m], NULL, ucl[[m]], uac[[m]]))) / M
    luy <- -.mlk(y, NULL, ycl, yac)
    lxy <- -.mlk(y, NULL, .cl, .ac)

    ## mean predictions
    mcv <- matrix(.0, N, N); .r <- 10
    for(.i in seq.int(.r))
    {
        .ik <- vapply(1:M, function(m) mvn(1, ucl[[m]], 1), y) %>% findKernel(ykn)
        mcv <- mcv + cmb(.ik, TUY)[[1]]
    }
    mcv <- mcv / .r                     # mean covariance of y
    .cl <- chol(mcv)
    .ac <- chol2inv(.cl)
    .ay <- drop(.ac %*% y)
    mnl <- nlk(y, NULL, .cl, .ac)       # mean negative log likelihood
    yht <- (mcv - diag(TUY[1], N)) %*% .ay
    ylo <- y - .ay / diag(.ac)
    mse <- mean((y - yht)^2)
    loo <- mean((y - ylo)^2)
    cyh <- cor(y, yht)
    
    ## update context and return
    for(. in ls()) .cx[[.]] <- get(.)
    .cx
}

## markder of derivatives
mdv <- c('!@#$%^&*') %>% strsplit(NULL) %>% unlist

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
sg3 <- function(ctx, bsz=30, lr0=1e-3, ...)
{
    ## list of history
    hst <- list(stt=list(), par=list())

    ## dots
    dot <- list(...)
    max.itr <- dot$max.itr %||% 1e4
    wtm <- dot$wtm %||% 5e-2
    tol <- dot$tol %||% 1e-9
    mmt <- dot$mmt %||% 0

    ## initial parameters
    lrt <- lr0
    mmt.dvt <- lapply(ctx$par, `*`, 0)
    i <- 1
    ## .lp <- rep(1:3, e=5) %>% rep(1000)
    .lp <- rep(1:3, l=1e5)
    .sk <- c('mnl', 'mse', 'cyh', 'lxy', 'luy', 'lxu')

    ## contex of lowest loss
    lw1 <- within(ctx, lxy <- Inf)
    lw2 <- within(ctx, mnl <- Inf)
    
    time.start <- Sys.time()
    while(TRUE)
    {
        tm0 <- Sys.time()
        ctx <- up3(ctx, .dv=.lp[i], ...)
        sec <- (Sys.time() - tm0) %T>% {units(.) <- 'secs'} %>% as.numeric
        mmt.dvt <- mapply(function(p, n) mmt * p + (1 - mmt) * n, mmt.dvt, ctx$dvt)
        
        ## update statistics upon new epoch
        for(. in .sk) assign(., ctx[[.]])
        
        ## record each iteration
        hst$stt[[i]] <- df(i=i, sec=sec, ctx[.sk])
        hst$par[[i]] <- ctx$par

        if(ctx$lxy < lw1$lxy)           # lowest Lxy
            lw1 <- within(ctx, {itr <- i})
        if(ctx$mnl < lw2$mnl)           # lowest Mnl
            lw2 <- within(ctx, {itr <- i})

        ## update learning rate
        if(i %% 80 == 1)
        {
            ## header of tracks
            hdr <- c(
                eph=sprintf("%5s", 'eph'),     # epoch
                lxy=sprintf("%7s", 'lxy'),
                luy=sprintf("%7s", 'luy'),
                lxu=sprintf("%7s", 'lxu'),
                mse=sprintf('%7s', 'mse'),     # developing erro
                cyh=sprintf('%4s', 'cyh'),     # developing corr
                mnl=sprintf('%7s', 'mnl'),     # covariance corr
                dvt=sprintf('%9s', 'xdv'),     # maximum derivative
                phi=sprintf('%7s', 'PHI'),     # phi
                sec=sprintf('%4s', 'sec'))     # seconds used
            cat(paste(hdr, collapse=' '), '\n', sep='')
        }

        ## message tracks
        xdv <- max(abs(unlist(mmt.dvt)))
        PHI <- exp(ctx$par$tuy[1, ])
        msg <- c(
            eph=sprintf("%04X",  i-1),  # epoch
            lxy=sprintf('%7.2f', lxy),  # error
            luy=sprintf('%7.2f', luy),  # error
            lxu=sprintf('%7.2f', lxu),  # error
            mse=sprintf('%7.2f', mse),  # error
            cyh=sprintf('%4.2f', cyh),  # corr(y, yhat)
            mnl=sprintf('%7.2f', mnl),  # corr(ycv.hat, ycv.ref)
            dvt=sprintf('%9.1e', xdv),  # maximum derivative
            phi=sprintf('%7.3f', PHI),  # phi
            sec=sprintf('%4.1f', sec))  # seconds used
        cat(mdv[.lp[i]], paste(msg, collapse=' '), '\n', sep='')
        if(i > max.itr)
        {
            cat('SGD: reaching max iter:', max.itr, '\n')
            break
        }
        if((Sys.time() - time.start) %T>% {units(.) <- 'hours'} %>% as.numeric > wtm && .lp[i] != 1)
        {
            cat('SGD: reaching walltime:', wtm, 'hour(s)\n'); break
            break
        }

        ## greadient descent here:
        new.par <- mapply(function(p, d) p - lrt * d, ctx$par, mmt.dvt)
        .df <- sum(abs(unlist(new.par) - unlist(ctx$par)))
        if(.df < tol && .lp[i] != 2)
            break
        ctx$par <- new.par
        i <- i + 1
    }

    ## rearrange history
    hst$stt <- do.call(rbind, hst$stt) %>% as.tbl %>% mutate(sec=cumsum(sec) - first(sec))

    ## some more aggregated statistics
    ctx$mcv <- mky(ctx, 1000); ctx$mnl <- with(ctx, nlk(y, mcv))
    lw1$mcv <- mky(lw1, 1000); lw1$mnl <- with(lw1, nlk(y, mcv))
    lw2$mcv <- mky(lw2, 1000); lw2$mnl <- with(lw2, nlk(y, mcv))
    
    ## return the history and new parameters
    ctx <- within(ctx, {hst <- hst; lw1 <- lw1; lw2 <- lw2})
    ctx
}

pl3 <- function(ctx, out=NULL, xlim=NULL, ylim=NULL, smooth=TRUE, ...)
{
    for(.n in names(ctx)) assign(.n, ctx[[.n]])

    ## reference LMM
    cyh.lmm <- rpt %>% filter(mtd=='lmm', dat=='dvp', key=='cyh') %>% select(val) %>% unlist
    mse.lmm <- rpt %>% filter(mtd=='lmm', dat=='dvp', key=='mse') %>% select(val) %>% unlist
    nlk.lmm <- rpt %>% filter(mtd=='lmm', dat=='dvp', key=='nlk') %>% select(val) %>% unlist

    ## phi
    PHI.sim <- exp(par$sim[1, ])
    PHI.hst <- sapply(hst$par, function(r) exp(r$tuy[1, ]))
    
    ## history relative to LMM
    rel <- hst$stt %>%
        mutate(mnl=mnl/nlk.lmm, mse=mse/mse.lmm, cyh=cyh/cyh.lmm,
               lxy=lxy/nlk.lmm, luy=luy/nlk.lmm, lxu=lxu/nlk.lmm, PHI=PHI.hst/PHI.sim)
    if(smooth)
        rel <- rel %>% mutate(mnl=lowess(sec, mnl, .1)[[2]],
                              mse=lowess(sec, mse, .1)[[2]], cyh=lowess(sec, cyh, .1)[[2]])
    ## rel <- rel %>% mutate_at(vars(mse:lxu), log)

    ## plot relative history
    if(!is.null(out))
        png(out, width=1024, height=1024, pointsize=24)
    cl <- rainbow(6)
    lt <- c(2, 2, 2, 1, 1, 1)
    if(is.null(xlim))
        xlim <- with(rel, range(sec))
    if(is.null(ylim))
        ylim <- with(rel, c(0, max(1, mse, cyh, mnl, lxy, luy, lxu)))
    within(rel,
    {
        plot(0, 0, 'n', xlim, ylim)
        lines(sec, mse, col=cl[1], lty=lt[1], lwd=2)
        lines(sec, cyh, col=cl[2], lty=lt[2], lwd=2)
        lines(sec, PHI, col=cl[3], lty=lt[3], lwd=2)
        lines(sec, mnl, col=cl[4], lty=lt[4], lwd=2)
        lines(sec, luy, col=cl[5], lty=lt[5], lwd=2)
        lines(sec, lxu, col=cl[6], lty=lt[6], lwd=2)

        abline(1, 0)

        legend(x='topright', legend=c('mse', 'cyh', 'PHI', 'lxy', 'luy', 'lxu'), lwd=2, lty=lt, col=cl)
    })
    if(!is.null(out))
        dev.off()
}

plc <- function(ctx, out=NULL)
{
    mcv.lw1 <- ctx$lw1$mcv
    mcv.lw2 <- ctx$lw2$mcv

    ycv.lw1 <- ctx$lw1$ycv
    ycv.lw2 <- ctx$lw2$ycv
    ycv.ref <- ctx$fit$ycv.ref

    par(mfrow=c(2, 3))
    image(mcv.lw1)
    image(mcv.lw2)
    image(ycv.ref)
    image(ycv.lw1)
    image(ycv.lw2)
    image(ycv.ref)
}
