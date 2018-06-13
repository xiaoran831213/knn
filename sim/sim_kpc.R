library(MASS)
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
#' @param gno list contains gnomic map, subject id, and genomic matrix
#' in alternative allele dosage format;
#' @param N numeric draw this many samples for model developing;
#' @param P numeric draw this many features (i.e., SNPs);
#' @param H numeric draw this many samples for model eveluation;
#' @param frq numeric percentage of functional features (i.e., casual SNPs);
#' @param ycv character the covariance structure to generate noise free outcome;
#' l: linear kernel (cross-product); G: Gaussian; L: Laplacian
#' @param eps numeric size of white noise adds to the noise free outcome;
#' @param M numeric the number of latent features;
#' @param ukn character the basic kernel types, for the kernel network;
#' l: linear kernel (cross-product)
#' @param ikn character the inner kernel types, for the kernel network;
#' features, required to calcuate the loss and gradient;
#' @param bsz numeric batch size for mini-batch based training, when NULL or
#' greater than the sample size N, the whole data training is resumed.
main <- function(gno, N, P, H=N, frq=.1, lnk=I, eps=.1, oks=c(id, p1), yks=c(id, p1), ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    ejt <- dot$ejt %||% .1
    arg <- match.call() %>% tail(-2) %>% as.list # %>% as.data.frame
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)
    
    if(is.character(gno)) gno <- readRDS(gno)
    if(is.null(gno)) gno <- readRDS('data/p35_cmn.rds')

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    gmx <- get.gmx(gno$gmx, N, H, P)
    vcs <- c(eps, runif(length(oks) - 1, 0, 3))

    ## simulation generated variables
    sim <- within(list(),
    {
        dvp <- get.sim(gmx$dvp, vcs, frq, lnk, oks, ejt)
        evl <- get.sim(gmx$evl, vcs, frq, lnk, oks, ejt)
    })

    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    ykn.dvp <- krn(gmx$dvp, yks, ...)   # kernel for training
    ykn.evl <- krn(gmx$evl, yks, ...)   # kernel for testing
    
    ## initialize parameters
    ini.txy <- matrix(rnorm(length(yks), sd=.05), length(yks), 1)

    ## use R's optimizer
    rop.dvp <- rop.lmm(sim$dvp$rsp, ykn.dvp, ini.txy)
    rop.evl.rpt <- knl.prd(sim$evl$rsp, ykn.evl, rop.dvp$par)
    rpt <- cl(rpt, DF(mtd='rop', dat='dvp', rop.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='rop', dat='evl', rop.evl.rpt))

    ## use polynomial MINQUE, order 1
    mq1.dvp <- knl.mnq(sim$dvp$rsp, ykn.dvp[-1], order=1)
    mq1.evl.rpt <- knl.mnq.evl(sim$evl$rsp, ykn.evl[-1], mq1.dvp$par, order=1)
    rpt <- cl(rpt, DF(mtd='mq1', dat='dvp', mq1.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='mq1', dat='evl', mq1.evl.rpt))
    
    ## use polynomial MINQUE, order 1, batched
    mq2.dvp <- kpc.mnq(sim$dvp$rsp, ykn.dvp[-1], sim$evl$rsp, ykn.evl[-1], order=1, ...)
    mq2.evl.rpt <- knl.mnq.evl(sim$evl$rsp, ykn.evl[-1], mq2.dvp$par, order=1)
    rpt <- cl(rpt, DF(mtd='mq2', dat='dvp', mq2.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='mq2', dat='evl', mq2.evl.rpt))

    ## use polynomial MINQUE, order 2
    mq3.dvp <- knl.mnq(sim$dvp$rsp, ykn.dvp[-1], order=2)
    mq3.evl.rpt <- knl.mnq.evl(sim$evl$rsp, ykn.evl[-1], mq3.dvp$par, order=2)
    rpt <- cl(rpt, DF(mtd='mq3', dat='dvp', mq3.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='mq3', dat='evl', mq3.evl.rpt))

    ## use polynomial MINQUE, order 2, batched
    mq4.dvp <- kpc.mnq(sim$dvp$rsp, ykn.dvp[-1], sim$evl$rsp, ykn.evl[-1], order=2, ...)
    mq4.evl.rpt <- knl.mnq.evl(sim$evl$rsp, ykn.evl[-1], mq4.dvp$par, order=2)
    rpt <- cl(rpt, DF(mtd='mq4', dat='dvp', mq4.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='mq4', dat='evl', mq4.evl.rpt))

    ## use GCTA:
    gct.dvp <- gcta.reml(sim$dvp$rsp, ykn.dvp)
    gct.evl.rpt <- knl.prd(sim$evl$rsp, ykn.evl, gct.dvp$par, logged=FALSE)
    rpt <- cl(rpt, DF(mtd='gct', dat='dvp', gct.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='gct', dat='evl', gct.evl.rpt))

    ## oracle fit and null fit:
    rpt <- cl(rpt, DF(sim$dvp$rpt, dat='dvp'))
    rpt <- cl(rpt, DF(sim$evl$rpt, dat='evl'))
    
    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}
