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
main <- function(N, P, Q=1, R=1, frq=.2, lnk=SC, eps=1, oks=~PL1, ...)
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
    dat <- get.sim(dat, frq, lnk, eps, oks, svc=svc, ...)
    
    ## training
    dvp <- with(dat$dvp,
    {
        ret <- list()
        kn1 <- krn(gmx, ~LN1)
        kn2 <- krn(gmx, ~JL2)
        ## ret <- CL(ret, IZM=IZM(rsp, knl, ...))
        ## ret <- CL(ret, OZ1=OZM(rsp, kn1, ...))
        ## ret <- CL(ret, ON1=ONM(rsp, kn1, ...))
        ## ret <- CL(ret, OZ2=OZM(rsp, kn2, ...))
        ret <- CL(ret, ON2=ONM(rsp, kn2, ...))
        ret <- CL(ret, GCT=GCT(rsp, kn1))
        ## ret <- CL(ret, FUL=GCT(rsp, krn(gmx,  oks)))
        ret <- CL(ret, NUL=NUL(rsp, NULL))
        ## ret <- CL(ret, JMX=PMQ(rsp, krn(gmx, ~JX1)))
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
        kn1 <- krn(gmx, ~LN1)
        kn2 <- krn(gmx, ~JL2)
        ## ret <- CL(ret, IZM=vpd(rsp, knl, dvp$IZM$par))
        ## ret <- CL(ret, OZ1=vpd(rsp, kn1, dvp$OZ1$par))
        ## ret <- CL(ret, ON1=vpd(rsp, kn1, dvp$ON1$par))
        ## ret <- CL(ret, OZ2=vpd(rsp, kn2, dvp$OZ2$par))
        ret <- CL(ret, ON2=vpd(rsp, kn2, dvp$ON2$par))
        ret <- CL(ret, GCT=vpd(rsp, kn1, dvp$GCT$par))
        ## ret <- CL(ret, FUL=vpd(rsp, krn(gmx,  oks), dvp$FUL$par))
        ret <- CL(ret, NUL=vpd(rsp, NULL, dvp$NUL$par))
        ## ret <- CL(ret, JMX=vpd(rsp, krn(gmx, ~JX1), dvp$JMX$par))
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
    r <- main(N=1024, P=10000, Q=3, R=1, vcs=1, mdl=AA, frq=.01)
}
