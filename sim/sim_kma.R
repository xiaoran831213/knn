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
main <- function(gno, N, P, Q=5, frq=.1, lnk=I, eps=.2, oks=c(id, p1), yks=c(id, p1), ...)
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
    if(is.null(gno)) gno <- readRDS('data/p35_cmn.rds')

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    nvc <- length(oks)
    gms <- sample.gmx(gno$gmx, N, P, Q + 1)
    vcs <- c(eps, rchisq(nvc - 1, 2))
    evl <- get.sim(gms[[Q + 1]], vcs, frq, lnk, oks, 0)
    dvp <- lapply(gms[1:Q], get.sim, vcs, frq, lnk, oks)
    
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    
    ## use polynomial MINQUE, order 2, batched, meta-analysis
    mq1 <- lapply(1:Q, function(i)
    {
        knl <- krn(gms[[i]], yks)
        ret <- kpc.mnq(dvp[[i]]$rsp, knl[-1], order=2, ...)
    })
    
    ## parameter estimates from meta-analysis
    par <- sapply(mq1, `[[`, 'par')     # par of all cohorts
    ac2 <- 1/sapply(mq1, `[[`, 'se2')   # ac2 of all cohorts
    ssz <- sapply(gms[1:Q], nrow)

    mq1.par <- rowSums(par*ac2) / rowSums(ac2) # by SE2
    mq2.par <- rowSums(par*ssz) / sum(ssz)     # by SSZ
    
    ## complete training sample
    gmx <- do.call(rbind, gms[1:Q])
    ykn <- krn(gmx, yks, ...)
    rsp <- do.call(c, lapply(dvp, '[[', 'rsp'))
    mq3 <- kpc.mnq(rsp, ykn[-1], order=2, ...) # MINQUE
    gct <- gcta.reml(rsp, ykn[-1])             # GCTA

    ## ----------------------- Testing Errors ----------------------- ##
    ykn <- krn(evl$gmx, yks, ...)
    rsp <- evl$rsp

    ## kernel MINQUE, batched, meta-analysis
    mq1.rpt <- knl.mnq.evl(rsp, ykn[-1], mq1.par, order=2, ...)
    rpt <- cl(rpt, DF(mtd='mq1', dat='evl', mq1.rpt))

    mq2.rpt <- knl.mnq.evl(rsp, ykn[-1], mq2.par, order=2, ...)
    rpt <- cl(rpt, DF(mtd='mq2', dat='evl', mq2.rpt))

    mq3.rpt <- knl.mnq.evl(rsp, ykn[-1], mq3$par, order=2, ...)
    rpt <- cl(rpt, DF(mtd='mq3', dat='evl', mq3.rpt))

    gct.rpt <- knl.prd(rsp, ykn, gct$par, logged=FALSE)
    rpt <- cl(rpt, DF(mtd='gct', dat='evl', gct.rpt))

    ## NULL
    rpt <- cl(rpt, DF(mtd='nul', dat='evl', nul(rsp)))

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}
