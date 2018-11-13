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
source('sim/eps.R')                     # noise

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
main <- function(N, P, Q=1, R=1, frq=.05, lnk=NL, eps=.1, oks=~LN1, ...)
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

    ## batched MINQUE: core algorithm
    assign('BMQ', dot$bmq %||% ONM, .GlobalEnv)
    assign('UBZ', dot$ubz %||%  64, .GlobalEnv)
    
    ## ------------------------- data genration ------------------------- ##
    ## for each of the Q groups, choose N samples and P features -> Training
    dat <- lapply(get.rds('sim/dat'), readRDS)
    dat <- get.gmx(dat, N, P, Q, R)
    dat <- get.sim(dat, frq=frq, lnk=lnk, eps=eps, oks=oks, ...)

    ## training
    dvp <- with(dat$dvp,
    {
        ret <- list()
        kn2 <- krn(gmx, ~JL2)
        kn1 <- krn(gmx, ~LN1)
        ret <- CL(ret, GCT=GCT(rsp, kn1))
        ret <- CL(ret, BM0=BM0(rsp, kn2, ...))
        ret <- CL(ret, BM1=BM1(rsp, kn2, ...))
        ret <- CL(ret, BM2=BM2(rsp, kn2, ...))
        ret <- CL(ret, BM3=BM3(rsp, kn2, ...))
        ret <- CL(ret, BMQ=BMQ(rsp, kn2, ...))
        ret <- CL(ret, NUL=NUL(rsp, NULL))
        ret
    })
    
    ## bias assesment
    vcs <- dat$dvp$vcs
    par <- pars(dvp, c(eps=eps, vcs))
    bia <- list(bias(dvp, eps, vcs))
    
    ## testing
    evl <- with(dat$evl,
    {
        ret <- list()
        kn2 <- krn(gmx, ~JL2)
        kn1 <- krn(gmx, ~LN1)
        ret <- CL(ret, GCT=vpd(rsp, kn1, dvp$GCT$par))
        ret <- CL(ret, BM0=vpd(rsp, kn2, dvp$BM0$par))
        ret <- CL(ret, BM1=vpd(rsp, kn2, dvp$BM1$par))
        ret <- CL(ret, BM2=vpd(rsp, kn2, dvp$BM2$par))
        ret <- CL(ret, BM3=vpd(rsp, kn2, dvp$BM3$par))
        ret <- CL(ret, BMQ=vpd(rsp, kn2, dvp$BMQ$par))
        ret <- CL(ret, NUL=vpd(rsp, NULL, dvp$NUL$par))
        ret
    })
    
    ## ----------------------- generate reports ----------------------- ##
    rpt <- list()
    dvp <- dvp %$% 'rpt'
    dvp <- lapply(names(dvp), function(.) cbind(dat='dvp', mtd=., dvp[[.]]))
    evl <- lapply(names(evl), function(.) cbind(dat='evl', mtd=., evl[[.]]))
    rpt <- c(dvp, evl, bia)
    
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
    r <- main(N=512, P=10000, Q=2, R=1, efn=EGS, eps=2.0, vcs=2, frq=.2, pss=0, bmq=ONM)
    r <- main(N=512, P=10000, Q=4, R=2, efn=EST, eps=1.5, vcs=3, frq=.2, oks=~PL, lnk=O2, ubz=64)
}
