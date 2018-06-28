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
main <- function(gno, N, P, H=N, frq=.1, lnk=I, eps=.2, oks=c(id, p1), yks=c(id, p1), ...)
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
    gmx <- sample.gmx(gno$gmx, N, P, Q=1, H=H) 
    vcs <- c(eps, rchisq(nvc - 1, 1))

    sim <- get.sim(gmx, vcs, frq, lnk, oks, ejt=c(rep(ejt, Q), 0.0))
    dvp <- sim[[1]]
    evl <- sim[[2]]
    gms <- lapply(dvp, `[[`, 'gmx')
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    
    ## use polynomial MINQUE, order 2, batched, meta-analysis
    mq1 <- {
        knl <- krn(dvp$gmx, yks)
        ret <- kpc.mnq(dvp$rsp, knl[-1], ...)}
    rop <- {
        knl <- krn(gms[[i]], yks)
        ret <- rop.lmm(dvp$rsp, knl)}
    
    ## ----------------------- Testing Errors ----------------------- ##
    ykn <- krn(evl$gmx, yks, ...)
    rsp <- evl$rsp

    mq4.rpt <- knl.mnq.evl(rsp, ykn[-1], mq4$par, order=2, ...)
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
