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
    gmx <- sample.gmx(gmx, N, P, Q, H) 
    vcs <- c(eps=eps, vc=rchisq(length(oks), 1))
    cvs <- c(eps=id, cv=oks)
    sim <- get.sim(gmx, vcs, frq, lnk, cvs, ejt=c(rep(ejt, Q), 0.0))
    dvp <- sim[1:Q]
    evl <- sim[[Q+1]]

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
        gmx <- do.call(rbind, lapply(dvp, `[[`, 'gmx'))
        rsp <- do.call(c, lapply(dvp, `[[`, 'rsp'))
        knl <- krn(gmx, yks)
        rpt <- list()

        ## average non meta analysis, 
        rpt <- cl(rpt, sapply(1:Q, function(i)
        {
            with(knl.prd(rsp, knl, sep$mnq[[i]]$fit$par, logged=FALSE, ...),
            {
                names(val) <- key; val
            })
        }) %>% rowMeans %>% {DF(mtd='avg.mnq', key=names(.), val=.)})
        rpt <- cl(rpt, sapply(1:Q, function(i)
        {
            with(knl.prd(rsp, knl, sep$rop[[i]]$fit$par, ...),
            {
                names(val) <- key; val
            })
        }) %>% rowMeans %>% {DF(mtd='avg.rop', key=names(.), val=.)})

        ## fit using all training samples
        mnq <- knl.mnq(rsp, knl, ...)
        rop <- rop.lmm(rsp, knl, ...)

        ## reports
        rpt <- cl(rpt, DF(mtd='whl.mnq', mnq$rpt))
        rpt <- cl(rpt, DF(mtd='whl.rop', rop$rpt))
        rpt <- cbind(dat='dvp', do.call(rbind, rpt))
        list(rpt=rpt, mnq=mnq, rop=rop)
    })

    ## ----------------------- Testing Errors ----------------------- ##
    evl <- with(evl,
    {
        knl <- krn(gmx, yks, ...)
        rpt <- list(rpt)
        ## meta analysis
        mnq <- agg.mat(sep$mnq)$a[, 'nlk']
        rop <- agg.mat(sep$rop)$a[, 'nlk']
        rpt <- cl(rpt, DF(mtd='mat.mnq', knl.prd(rsp, knl, mnq, logged=FALSE, ...)))
        rpt <- cl(rpt, DF(mtd='mat.rop', knl.prd(rsp, knl, rop, ...)))
        
        ## average non meta analysis
        rpt <- cl(rpt, sapply(1:Q, function(i)
        {
            with(knl.prd(rsp, knl, sep$mnq[[i]]$fit$par, logged=FALSE, ...),
            {
                names(val) <- key; val
            })
        }) %>% rowMeans %>% {DF(mtd='avg.mnq', key=names(.), val=.)})
        rpt <- cl(rpt, sapply(1:Q, function(i)
        {
            with(knl.prd(rsp, knl, sep$rop[[i]]$fit$par, ...),
            {
                names(val) <- key; val
            })
        }) %>% rowMeans %>% {DF(mtd='avg.rop', key=names(.), val=.)})

        ## whole data model (golden standard)
        mnq <- dvp$mnq$par
        rop <- dvp$rop$par
        rpt <- cl(rpt, DF(mtd='whl.mnq', knl.prd(rsp, knl, mnq, logged=FALSE, ...)))
        rpt <- cl(rpt, DF(mtd='whl.rop', knl.prd(rsp, knl, rop, ...)))
        rpt=cbind(dat='evl', do.call(rbind, rpt))
        list(rpt=rpt)
    })
    ## mq4.rpt <- knl.mnq.evl(rsp, ykn[-1], mq4$par, order=2, ...)
    ## rpt <- cl(rpt, DF(mtd='mq4', dat='evl', mq4.rpt))
    
    ## gct.rpt <- knl.prd(rsp, ykn, gct$par, logged=FALSE)
    ## rpt <- cl(rpt, DF(mtd='gct', dat='evl', gct.rpt))

    rpt <- list(dvp=dvp$rpt, evl=evl$rpt)

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}

#' aggregate cohorts
agg.mat <- function(rpt, ...)
{
    ## leave one cohort out validation
    w <- list(length(rpt))
    p <- list(length(rpt))
    for(i in seq_along(rpt))
    {
        org <- rpt[[i]]                 # the i.th cohort
        pol <- rpt[i]                   # exclude the i th.
        fit <- org[['fit']]         
        mtv <- lapply(pol, function(dat) # mata validation
        {
            with(knl.prd(dat$rsp, dat$knl, fit$par, ...),
            {
                names(val) <- key
                val
            })
        })
        mtv <- as.data.frame(do.call(rbind, mtv))
        mtv <- with(mtv,
        {
            c(mse=mean(1/mse),
              nlk=mean(1/nlk),
              cyh=mean(cyh),
              ssz=NROW(org$rsp))
        })
        w[[i]] <- mtv
        p[[i]] <- fit[['par']]
    }

    ## row: parameter, col: cohort
    p <- do.call(cbind, p)

    ## row: cohort, col: error-measurement
    w <- as.data.frame(do.call(rbind, w))
    w <- sapply(w, function(.) . / sum(.), simplify=FALSE)
    w <- do.call(cbind, w)

    ## weighted aggregation
    a <- p %*% w
    colnames(a) <- colnames(w)

    list(p=p, w=w, a=a)
}


test <- function()
{
    r <- main(N=200, P=2000, Q=5, H=1000, frq=.1, eps=.1, oks=c(p1), yks=c(p1), ejt=0.2)
}
