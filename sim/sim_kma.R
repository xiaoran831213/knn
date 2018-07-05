library(magrittr)
library(Matrix)
source('R/hlp.R')
source('R/kpl.R')
source('R/gct/grm.R')
source('R/gct/gct.R')
source('R/utl.R')
source('R/lmm.R')
source("R/mnq.R")
source("R/kmq.R")
source("sim/sim_kpl.R")
source('sim/utl.R')

library(devtools)
devtools::load_all()

#' simulation of kernel deep neural network;
#' @param gno gnomic map, subject id, and genomic matrix in dosage format;
#' @param N numeric draw this many samples for model developing;
#' @param P numeric draw this many features (i.e., SNPs);
#' @param H numeric draw this many samples for model eveluation;
#' @param frq numeric percentage of functional features (i.e., casual SNPs);
#' @param eps numeric size of white noise adds to the noise free outcome;
#' @param bsz numeric batch size for mini-batch based training, when NULL or
main <- function(N, P, Q=5, H=N, frq=.1, lnk=I, eps=.2, oks=c(p1), yks=c(p1), ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    ejt <- dot$ejt %||% .2
    arg <- match.call() %>% tail(-1) %>% as.list # %>% as.data.frame
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    gmx <- readRDS('data/p35_c05.rds')$gmx
    gmx <- sample.gmx(gmx, N, P, Q)
    vcs <- c(eps=eps, vc=rchisq(length(oks), 1))
    cvs <- c(eps=id, cv=oks)
    ## sim <- get.sim(gmx, vcs, frq, lnk, cvs, ejt=c(rep(ejt, Q), 0.0))
    dvp <- get.sim(gmx$dvp, vcs, frq, lnk, cvs, ejt=ejt)
    evl <- get.sim(gmx$evl, vcs, frq, lnk, cvs, ejt=0.0)

    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()

    ## fit each cohort separately
    mnq <- lapply(1:Q, function(i)
    {
        within(dvp[[i]], {knl <- krn(gmx, yks); fit <- knl.mnq(rsp, knl, ...)})
    })
    rop <- lapply(1:Q, function(i)
    {
        within(dvp[[i]], {knl <- krn(gmx, yks); fit <- rop.lmm(rsp, knl, ...)})
    })
    sep <- list(mnq=lapply(mnq, `[`, c('rsp', 'fit', 'knl')),
                rop=lapply(rop, `[`, c('rsp', 'fit', 'knl')))

    dvp <- with(list(),
    {
        gmx <- do.call(rbind, lapply(dvp, `[[`, 'gmx')) # combined G
        rsp <- unlist(lapply(dvp, `[[`, 'rsp'))         # combined Y
        knl <- krn(gmx, yks)                            # combined K
        rpt <- list()
        
        ## average non meta analysis, 
        mnq <- sapply(sep$mnq, function(.) knl.prd(rsp, knl, .$fit$par, ln=0, rt=0)) %>%
            rowMeans %>% {DF(mtd='mnq.avg', key=names(.), val=.)}
        rop <- sapply(sep$rop, function(.) knl.prd(rsp, knl, .$fit$par, ln=0, rt=0)) %>%
            rowMeans %>% {DF(mtd='rop.avg', key=names(.), val=.)}
        rpt <- CL(rpt, mnq, rop)
        
        ## fit using all training samples
        mnq <- knl.mnq(rsp, knl, ...)
        rop <- rop.lmm(rsp, knl, ...)
        rpt <- CL(rpt, DF(mtd='whl.mnq', mnq$rpt), DF(mtd='whl.rop', rop$rpt))

        rpt <- cbind(dat='dvp', do.call(rbind, rpt))
        list(rpt=rpt, mnq=mnq, rop=rop)
    })

    ## ----------------------- Testing Errors ----------------------- ##
    evl <- with(list(),
    {
        gmx <- do.call(rbind, lapply(evl, `[[`, 'gmx')) # combined G
        rsp <- unlist(lapply(evl, `[[`, 'rsp'))         # combined Y
        knl <- krn(gmx, yks, ...)
        rpt <- list()
        ## performance: meta analysis
        mnq <- agg.mat(sep$mnq)$a[, 'nlk']
        rpt <- cl(rpt, DF(mtd='nlk.mnq', knl.prd(rsp, knl, mnq, ln=0, ...)))
        mnq <- agg.mat(sep$mnq)$a[, 'mse']
        rpt <- cl(rpt, DF(mtd='mse.mnq', knl.prd(rsp, knl, mnq, ln=0, ...)))
        rop <- agg.mat(sep$rop)$a[, 'nlk']
        rpt <- cl(rpt, DF(mtd='nlk.rop', knl.prd(rsp, knl, rop, ln=0, ...)))

        ## performance: average of non meta analysis
        mnq <- sapply(1:Q, function(i)
        {
            knl.prd(rsp, knl, sep$mnq[[i]]$fit$par, ln=0, rt=0, ...)
        })
        med <- apply(t(mnq), 2L, median)
        med <- DF(mtd='med.mnq', key=names(med), val=med)
        avg <- apply(t(mnq), 2L, mean)
        avg <- DF(mtd='avg.mnq', key=names(avg), val=avg)
        rpt <- CL(rpt, med, avg)

        rop <- sapply(1:Q, function(i)
        {
            knl.prd(rsp, knl, sep$rop[[i]]$fit$par, ln=0, rt=0, ...)
        })
        med <- apply(t(rop), 2L, median)
        med <- DF(mtd='med.rop', key=names(med), val=med)
        avg <- apply(t(rop), 2L, mean)
        avg <- DF(mtd='avg.rop', key=names(avg), val=avg)
        rpt <- CL(rpt, med, avg)
        
        ## performance: model of whole training data
        mnq <- dvp$mnq$par
        rop <- dvp$rop$par
        rpt <- CL(rpt, DF(mtd='whl.mnq', knl.prd(rsp, knl, mnq, ln=0, ...)))
        rpt <- CL(rpt, DF(mtd='whl.rop', knl.prd(rsp, knl, rop, ln=0, ...)))

        ## report
        rpt=cbind(dat='evl', do.call(rbind, rpt))
        list(rpt=rpt)
    })
    rpt <- list(dvp=dvp$rpt, evl=evl$rpt)

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}

#' leave one out assessment
loo.mat <- function(rpt, ...)
{
    ret <- lapply(seq_along(rpt), function(i)
    {
        mtv <- sapply(seq_along(rpt)[-i], function(j) # mata validation
        {
            . <- rpt[[j]]
            . <- knl.prd(.$rsp, .$knl, .$fit$par, rt=0, ...)
            c(j=j, .)
        })
        DF(i=i, t(mtv))
    })
    do.call(rbind, ret)
}


#' aggregate cohorts
agg.mat <- function(rpt, mat=0, ...)
{
    ## report of leave one (cohort) out
    ## row: cohort, col: errors
    r <- loo.mat(rpt, ...)
    f <- r$i
    r <- subset(r, se=-c(i, j))

    ## row: parameter, col: cohort
    p <- sapply(rpt, function(.) .$fit$par)

    ## row: cohort, col: weights
    if(mat==0)
    {
        w <- within(r,
        {
            nlk <- 1/nlk
            mse <- 1/mse
            cyh <- 1/(1-cyh)
        })
        w <- by(w, f, function(.)
        {
            sapply(., mean)
        })
        w <- do.call(rbind, w)
    }
    else
    {
        w <- by(r, f, function(.)
        {
            sapply(., mean)
        })
        w <- as.data.frame(do.call(rbind, w))
        w <- within(w,
        {
            nlk <- 1/nlk
            mse <- 1/mse
            cyh <- 1/(1-cyh)
        })
        w <- as.matrix(w)
    }
    w <- sweep(w, 2L, colSums(w), '/')
    ## weighted aggregation
    a <- p %*% w

    list(p=p, w=w, a=a)
}

test <- function()
{
    r <- main(N=200, P=2000, Q=5, H=1000, frq=.1, eps=.1, oks=c(p1), yks=c(p1), ejt=0.2)
}
