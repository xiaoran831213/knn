library(magrittr)
library(Matrix)

source('R/hlp.R')                       # helpers
source('R/kpl.R')                       # kernel players
source('R/utl.R')                       # utilities
source('R/vcm.R')                       # variance component models (VCM)
source('R/mnq.R')                       # MINQUE
source('R/msg.R')                       # message board
source('R/bat.R')                       # batched VCM trainer
source('R/agg.R')                       # model aggregation
source('R/bat.R')
source('sim/grm.R')                     # genomic relatedness matrix (plink, GCTA)
source('sim/gct.R')                     # GCTA wrapper
source('sim/utl.R')                    # simulation utilites
source('sim/gsm.R')                     # simulation utilites

library(devtools)                       # enable the C++ functions
devtools::load_all()

#' simulation of kernel deep neural network;
#' @param N size of population groups
#' @param P number of variants
#' @param Q number of groups that constitute training data
#' @param R number of groups that constitute testing data
#' @param frq fraction of functional variants
#' @param eps size of white noise
#' 
#' @param oks true kernels to generate y - the responses;
#' @param svc variance of true variance components;
#' @param lnk link function to transform the generated response;
#'
#' @param yks used kernels for modeling;
#' 
#' @param bsz batch size for batched training
#'
#' see "sim/utl.R" to understand {oks}, {lnk}, and {yks}.
main <- function(N, P, Q=1, R=0, frq=.05, lnk=I, eps=.1, oks=p1, yks=p1, ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    het <- dot$het %||% 0.0
    svc <- dot$svc %||% 1.0
    vcs <- dot$vcs %||% get.vcs(length(oks), 'r', svc)
    arg <- match.call() %>% tail(-1) %>% as.list
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)
    gds <- dot$gds %||% 'ukb'
    if(gds=='ukb')
        gds <- get.rds('sim/dat')
    else
        gds <- 'data/1kg_c05.rds'
    
    ## ------------------------- data genration ------------------------- ##
    ## for each of the Q groups, choose N samples and P features -> Training
    gmx <- get.gmx(readRDS(gds), N, P, Q, R)
    dat <- with(gmx, c(dvp, evl))
    dat <- get.sim(dat, frq, lnk, eps, vcs, oks, ...)
    dvp <- dat[+(1:Q)]
    gmx <- do.call(rbind, dvp %$% 'gmx') # stacked genomics
    rsp <- unlist(dvp %$% 'rsp')         # stacked response
    knl <- krn(gmx, yks)                # kernels from the entire data

    ## ------------------------- model fit ------------------------- ##
    dvp <- within(list(),
    {
        gct <- gcta.reml(rsp, knl)
        mle <- rop.vcm(rsp, knl)
        bml <- GBT(rop.vcm, rsp, knl, ...)
        mnq <- knl.mnq(rsp, knl, prd=FALSE, psd=FALSE)
        bmq <- GBT(knl.mnq, rsp, knl, prd=FALSE, psd=FALSE, ...)
    })

    ## bias assesment
    par <- dvp %$% 'par'                # estimates
    ref <- c(eps, vcs)                  # reference
    len <- max(length(ref), do.call(pmax, par %$% length))
    key <- paste0('bia.', c('eps', paste0('vc', 1:(len-1))))

    ## zero padding, and biases
    par <- lapply(par, function(.) c(., rep(0, len - length(.))))
    ref <- c(ref, rep(0, len - length(ref)))
    bia <- lapply(names(par), function(.) DF(mtd=., key=key, val=par[[.]] - ref))
    bia <- DF(dat='dvp', do.call(rbind, bia))
    
    ## ----------------------- generate reports ----------------------- ##
    rpt <- list()
    rpt <- CL(rpt, bia)

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    print(c(par, list(ref=ref)))
    ret
}

test <- function()
{
    r <- main(N=500, P=1000, Q=2, R=0, frq=.1, eps=1, oks=p2, yks=p2, svc=1)
    subset(r, dat=='evl' & key=='nlk' | dat=='dvp' & key=='rtm')
}
