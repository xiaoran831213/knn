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
source('R/bat.R')                       # batched trainer
source('sim/grm.R')                     # genomic relatedness matrix (plink, GCTA)
source('sim/gct.R')                     # GCTA wrapper
source('sim/utl.R')                     # simulation utilites
source('sim/mdl.R')                     # models
source('sim/lnk.R')                     # link functions
source('sim/gsm.R')                     # genomic simulator
source('sim/mtd.R')                     # methods

library(devtools)                       # enable the C++ functions
devtools::load_all()

## methodology
mkm <- function(mtd=knl.mnq, mdl=~LN1, ...)
{
    list(mtd=mtd, mdl=mdl, arg=list(...))
}
MTD <- within(list(),
{

})

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
main <- function(N, P, Q=1, R=1, frq=.05, lnk=I, eps=.1, rks=~LN1, oks=~p1+ga, yks=~p1, ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    het <- dot$het %||% 0.0
    svc <- dot$svc %||% 1.0
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
    dat <- lapply(gds, readRDS)
    dat <- get.gmx(dat, N, P, Q, R)
    dat <- get.sim(dat, frq, lnk, eps, oks, yks, svc=svc, ...)

    ## training
    dvp <- with(dat$dvp,
    {
        ret <- list()
        ret <- CL(ret, gct=GCT(rsp, krn(gmx, rks)))
        ret <- CL(ret, ful=GCT(rsp, krn(gmx, oks)))
        ret <- CL(ret, pmq=PMQ(rsp, knl))
        ret <- CL(ret, bm2=BM2(rsp, knl, ...))
        ret <- CL(ret, bm3=BM3(rsp, knl, ...))
        ret <- CL(ret, bm4=BM4(rsp, knl, ...))
        ret <- CL(ret, nul=NUL(rsp, knl))
        ret
    })
    
    ## bias assesment
    vcs <- dat$dvp$vcs
    par <- CL(dvp %$% 'par', ref=c(eps, vcs))
    bia <- bias(dvp, eps, vcs)

    ## testing
    evl <- with(dat$evl,
    {
        ret <- list()
        ret <- CL(ret, gct=vpd(rsp, krn(gmx, rks), dvp$gct$par))
        ret <- CL(ret, ful=vpd(rsp, krn(gmx, oks), dvp$ful$par))
        ret <- CL(ret, pmq=vpd(rsp, knl, dvp$pmq$par))
        ret <- CL(ret, bm2=vpd(rsp, knl, dvp$bm2$par))
        ret <- CL(ret, bm3=vpd(rsp, knl, dvp$bm3$par))
        ret <- CL(ret, bm4=vpd(rsp, knl, dvp$bm4$par))
        ret <- CL(ret, nul=vpd(rsp, NULL, dvp$nul$par))
        ret
    })
    
    ## ----------------------- generate reports ----------------------- ##
    rpt <- list()
    dvp <- dvp %$% 'rpt'
    dvp <- lapply(names(dvp), function(.) cbind(dat='dvp', mtd=., dvp[[.]]))
    evl <- lapply(names(evl), function(.) cbind(dat='evl', mtd=., evl[[.]]))
    rpt <- c(dvp, evl, list(bia))
    
    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 5L))
    ret <- cbind(arg, rpt)

    print(list(par=par))
    print(list(time=subset(ret, dat=='dvp' & key=='rtm')))
    ret
}

test <- function()
{
    u <- main(N=500, P=10000, Q=2, R=1, frq=.1, eps=.5, svc=2, oks=~LN2, yks=~LN2)
    subset(r, dat=='evl' & key=='nlk' | dat=='dvp' & key=='rtm')
}
