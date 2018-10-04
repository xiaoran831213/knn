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

library(devtools)                       # enable the C++ functions
devtools::load_all()

## methodology
## {response model, data model, data}
##  

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
        ret <- CL(ret, gct=gct.rml(rsp, krn(gmx, rks)))
        ## ret <- CL(ret, mle=rop.vcm(rsp, knl))
        ret <- CL(ret, sat=gct.rml(rsp, krn(gmx, oks)))
        ret <- CL(ret, mn0=knl.mnq(rsp, knl, psd=FALSE))
        ret <- CL(ret, mn1=knl.mnq(rsp, knl, psd=TRUE))
        ret <- CL(ret, nul=nul.vcm(rsp))
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
        ## ret <- CL(ret, mle=vpd(rsp, knl, dvp$mle$par))
        ret <- CL(ret, sat=vpd(rsp, krn(gmx, oks), dvp$sat$par))
        ret <- CL(ret, mn0=vpd(rsp, knl, dvp$mn0$par))
        ret <- CL(ret, mn1=vpd(rsp, knl, dvp$mn1$par))
        ret <- CL(ret, nul=nul.vcm(rsp, dvp$nul$par)$rpt)
        ## ret <- CL(ret, fun=vpd(rsp, fnl, c(eps, vcs)))
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
