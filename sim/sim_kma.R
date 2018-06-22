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
    gms <- sample.gmx(gno$gmx, N, P, Q)
    vcs <- sample.vcs(eps, nvc, Q)

    sim.dvp <- mapply(function(g, v)
    {
        get.sim(g, v, frq, lnk, oks, ejt)
    },
    gms, vcs, SIMPLIFY=FALSE)
    names(sim) <- sprintf('S%d', seq_along(sim))
    
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    
    ## use polynomial MINQUE, order 2, batched, meta-analysis
    mq1 <- mapply(function(C, K)
    {
        kpc.mnq(C$rsp, K[-1], order=2, ...)
    },
    sim, lapply(gms, krn, yks, ...), SIMPLIFY=FALSE)

    ## parameter estimates from meta-analysis
    par <- sapply(mq1, `[[`, 'par')     # par of all cohorts
    ac2 <- 1/sapply(mq1, `[[`, 'se2')   # ac2 of all cohorts
    mq1.par <- colSums(par*ac2) / colSums(ac2) # by SE2
    mq2.par <- colMeans(par)                   # by SSZ

    ## complete sample
    gmx <- do.call(rbind, gms)
    ykn <- krn(gmx, yks, ...)
    rsp <- do.call(c, lapply(sim, '[[', 'rsp'))

    ## kernel MINQUE, batched, complete
    mq1.rpt <- knl.mnq.evl(rsp, ykn[-1], mq1.par, order=2, ...)
    rpt <- cl(rpt, DF(mtd='mq1', dat='sep', mq1.rpt))

    mq2.rpt <- knl.mnq.evl(rsp, ykn[-1], mq2.par, order=2, ...)
    rpt <- cl(rpt, DF(mtd='mq2', dat='sep', mq2.rpt))

    ## kernel MINQUE, batched, meta-analysis
    mq3 <- kpc.mnq(rsp, ykn[-1], order=2, ...)
    rpt <- cl(rpt, DF(mtd='mq3', dat='cmb', mq3$rpt))
    
    ## GCTA
    gct <- gcta.reml(rsp, ykn)
    rpt <- cl(rpt, DF(mtd='gct', dat='cmb', gct$rpt))

    ## NULL
    rpt <- cl(rpt, DF(mtd='nul', dat='cmb', nul(rsp)))

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}
