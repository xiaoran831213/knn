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
main <- function(gno, N, P, Q=5, H=N, mat='s2',
                 frq=.1, lnk=I, eps=.2, oks=c(id, p1), yks=c(id, p1), ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    ejt <- dot$ejt %||% .3
    arg <- match.call() %>% tail(-2) %>% as.list # %>% as.data.frame
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)

    if(is.character(gno)) gno <- readRDS(gno)
    if(is.null(gno)) gno <- readRDS('data/p35_c05.rds')

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    nvc <- length(oks)
    gmx <- sample.gmx(gno$gmx, N, P, Q, H) 
    vcs <- c(eps, rchisq(nvc - 1, 1))
    ## vcs <- c(eps, rep(2.0, nvc - 1))
    sim <- get.sim(gmx, vcs, frq, lnk, oks, ejt=c(rep(ejt, Q), 0.0))
    dvp <- sim[1:Q]
    evl <- sim[[Q+1]]
    gms <- lapply(dvp, `[[`, 'gmx')
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    
    ## use polynomial MINQUE, order 2, batched, meta-analysis
    mq1 <- lapply(1:Q, function(i)
    {
        knl <- krn(gms[[i]], yks)
        ret <- kpc.mnq(dvp[[i]]$rsp, knl[-1], order=2, ...)
    })
    rop <- lapply(1:Q, function(i)
    {
        knl <- krn(gms[[i]], yks)
        ret <- rop.lmm(dvp[[i]]$rsp, knl)
    })
    
    ## complete training sample
    gmx <- do.call(rbind, gms[1:Q])
    ykn <- krn(gmx, yks, ...)
    rsp <- do.call(c, lapply(dvp, '[[', 'rsp'))
    ## mq3 <- kpc.mnq(rsp, ykn[-1], order=2, ...) # MINQUE mini-batched
    mq4 <- knl.mnq(rsp, ykn[-1], order=2) # MINQUE whole sample
    ## gct <- gcta.reml(rsp, ykn[-1])             # GCTA

    ## ----------------------- Testing Errors ----------------------- ##
    ykn <- krn(evl$gmx, yks, ...)
    rsp <- evl$rsp

    ## kernel MINQUE, batched, meta-analysis
    mq1.rpt <- lapply(1:Q, function(i)
    {
        . <- mq1[1:i]
        p <- sapply(., `[[`, 'par')     # par of all cohorts
        if(mat == 's2')
        {
            a <- 1/sapply(., `[[`, 'se2') # ac2 of all cohorts
            p <- rowSums(p * a) / rowSums(a) # by SE2
        }
        if(mat == 'sz')
        {
            a <- sapply(gms[1:i], nrow)
            p <- rowSums(sweep(p, 2L, a, '*')) / sum(a) # by SSZ
        }
        if(mat == 'me')
        {
            a <- sapply(., function(x) unlist(subset(x$rpt, key=='mse', 'val')))
            a <- 1/sqrt(a)
            p <- rowSums(sweep(p, 2L, a, '*')) / sum(a) # by MSE
        }
        if(mat == 'lk')
        {
            a <- 1/sapply(., function(x) unlist(subset(x$rpt, key=='nlk', 'val')))
            p <- rowSums(sweep(p, 2L, a, '*')) / sum(a) # by nlk
        }
        print(a)

        m <- sprintf('m%02d', i)
        DF(mtd=m, dat='evl', knl.mnq.evl(rsp, ykn[-1], p, order=2, ...))
    })
    rpt <- cl(rpt, do.call(rbind, mq1.rpt))

    ## mq4.rpt <- knl.mnq.evl(rsp, ykn[-1], mq4$par, order=2, ...)
    ## rpt <- cl(rpt, DF(mtd='mq4', dat='evl', mq4.rpt))
    
    ## gct.rpt <- knl.prd(rsp, ykn, gct$par, logged=FALSE)
    ## rpt <- cl(rpt, DF(mtd='gct', dat='evl', gct.rpt))

    rop.rpt <- lapply(1:Q, function(i)
    {
        p <- rowMeans(sapply(rop[1:i], `[[`, 'par'))
        m <- sprintf('r%02d', i)
        DF(mtd=m, knl.prd(rsp, ykn, p))
    })

    rop.rpt <- do.call(rbind, rop.rpt)
    rpt <- cl(rpt, DF(dat='evl', rop.rpt))
    
    ## NULL
    rpt <- cl(rpt, DF(mtd='nul', dat='evl', nul(rsp)))
    ## GOLD
    rpt <- cl(rpt, DF(dat='evl', evl$rpt))

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}

test <- function()
{
    r <- main(NULL, N=500, P=2500, Q=2, H=1500, frq=.1,  eps=.1, oks=c(id, p1), yks=c(id, p1), ejt=0.1, bsz=100, wep=2)
}
