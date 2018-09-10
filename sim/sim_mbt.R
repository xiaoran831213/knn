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
source('sim/grm.R')                     # genomic relatedness matrix (plink, GCTA)
source('sim/gct.R')                     # GCTA wrapper
source('sim/utl.R')                     # simulation utilites
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
main <- function(N, P, Q=1, R=1, frq=.05, lnk=I, eps=.1, oks=p1, yks=p1, ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    het <- dot$het %||% .0
    svc <- dot$svc %||% 1
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
    dat <- get.sim(dat, frq, lnk, eps, oks, het=c(rep(het, Q), rep(het, R)), ...)
    dvp <- dat[+(1:Q)]                  # training
    evl <- dat[-(1:Q)]                  # testing

    ## ------------------------- model building ------------------------- ##
    dvp <- within(list(),
    {
        gmx <- do.call(rbind, EL2(dvp, 'gmx')) # stacked genomics
        rsp <- unlist(EL2(dvp, 'rsp'))         # stacked response
        knl <- krn(gmx, yks)                   # kernels from the entire data
        gct <- gcta.reml(rsp, krn(gmx, p1))    # whole sample GCTA-REML
        mnq <- knl.mnq(rsp, knl, cpp=TRUE)     # whole sample MINQUE

        ## batched MINQUE
        kmq <- GBT(knl.mnq, rsp, knl, cpp=TRUE, ...)
        kmq.par <- kmq$par
        kmq.rtm <- kmq$rtm
        kmq <- do.call(rbind, lapply(kmq.par, vpd, y=rsp, K=knl))
        kmq.rpt <- DF(mtd=paste0('mnq.', sub("[.][^.]*$", "", rownames(kmq))), kmq)
        kmq.rpt <- rbind(kmq.rpt, DF(mtd='mnq.bat', key='rtm', val=kmq.rtm))
    })

    ## ------------------------- model evaluation ------------------------- ##
    evl <- within(list(),
    {
        gmx <- do.call(rbind, EL2(evl, 'gmx')) # stacked genomics
        rsp <- unlist(EL2(evl, 'rsp'))         # stacked response
        knl <- krn(gmx, yks)                   # kernels from the entire data

        gct <- DF(mtd='gct', vpd(rsp, krn(gmx, p1), dvp$gct$par)) # evaluate GCTA-REML
        mnq <- DF(mtd='mnq.whl', vpd(rsp, knl, dvp$mnq$par))      # evaluate MINQUE

        ## batched MINQUE
        kmq <- do.call(rbind, lapply(dvp$kmq.par, vpd, y=rsp, K=knl))
        kmq <- DF(mtd=paste0('mnq.', sub("[.][^.]*$", "", rownames(kmq))), kmq)
    })
    
    ## --------------------------- make reports --------------------------- ##
    rpt <- list()
    rpt <- CL(rpt, DF(dat='dvp', mtd='mnq.whl', dvp$mnq$rpt))
    ## rpt <- CL(rpt, DF(dat='dvp', mtd='mle.whl', dvp$mle$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='gct', dvp$gct$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='nul', nul(dvp$rsp)))

    rpt <- CL(rpt, DF(dat='evl', evl$mnq))
    ## rpt <- CL(rpt, DF(dat='evl', evl$mle))
    rpt <- CL(rpt, DF(dat='evl', evl$gct))
    rpt <- CL(rpt, DF(dat='evl', mtd='nul', nul(evl$rsp)))

    rpt <- CL(rpt, DF(dat='dvp', dvp$kmq.rpt)) # batched
    rpt <- CL(rpt, DF(dat='evl', evl$kmq))

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    print(list(kmq=dvp$kmq$par, mnq=dvp$mnq$par, gct=dvp$gct$par))
    invisible(ret)
}

test <- function()
{
    r <- main(N=2000, P=4000, Q=1, R=1, eps=.01, yks=p2, svc=5, mdl=ax, bsz=200)
    subset(r, dat == 'evl' & key =='nlk')
}
