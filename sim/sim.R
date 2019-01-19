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

## library(devtools)                       # enable the C++ functions
## devtools::load_all()

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
    efn <- dot$efn %||% EGS             # epsilon function for noise
    arg <- as.list(match.call()[-1])
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)
    gds <- dot$gds %||% 'ukb'
    gds <- if(gds=='ukb') get.rds('sim/dat') else 'data/1kg_c05.rds'
    
    ## ------------------------- data genration ------------------------- ##
    ## for each of the Q groups, choose N samples and P features -> Training
    dat <- lapply(gds, readRDS)
    dat <- get.gmx(dat, N, P, Q, R)
    dat <- get.sim(dat, frq=frq, lnk=lnk, eps=eps, oks=oks, ...)

    std <- function(k) k / mean(diag(k))
    ## training
    dvp <- with(dat$dvp,
    {
        ret <- list()
        ## kn1 <- krn(gmx, ~LN1)
        ## kn2 <- krn(gmx, ~LN2)
        kn3 <- lapply(krn(gmx, ~LN3), std)
        kn1 <- kn3[1]
        kn2 <- kn3[1:2]
        ret <- CL(ret, NUL=NUL(rsp, ...))
        ret <- CL(ret, G1K=GCT(rsp, kn1))
        ## ret <- CL(ret, M1K=MQ3(rsp, kn1, ...))
        ret <- CL(ret, BM2=MQ3(rsp, kn2, ...))
        ret <- CL(ret, BM3=MQ3(rsp, kn3, ...))
        ## ret <- CL(ret, UM2=MQ0(rsp, kn2, ...))
        ## ret <- CL(ret, UM3=MQ0(rsp, kn3, ...))
        ret
    })

    ## testing
    evl <- with(dat$evl,
    {
        ret <- list()
        ## kn1 <- krn(gmx, ~LN1)
        kn3 <- lapply(krn(gmx, ~LN3), std)
        kn1 <- kn3[1]
        kn2 <- kn3[1:2]
        ret <- CL(ret, NUL=vpd(dvp$NUL$par, rsp, rt=0))
        ret <- CL(ret, G1K=vpd(dvp$G1K$par, rsp, kn1, rt=0))
        ## ret <- CL(ret, M1K=vpd(dvp$M1K$par, rsp, kn1, X=x00, rt=0))
        ret <- CL(ret, BM2=vpd(dvp$BM2$par, rsp, kn2, rt=0))
        ret <- CL(ret, BM3=vpd(dvp$BM3$par, rsp, kn3, rt=0))
        ## ret <- CL(ret, UM2=vpd(dvp$UM2$par, rsp, kn2, rt=0))
        ## ret <- CL(ret, UM3=vpd(dvp$UM3$par, rsp, kn3, rt=0))
        ret
    })
    
    ## ----------------------- generate reports ----------------------- ##
    rtm <- unlist(dvp %$% 'rtm')
    par <- do.call(.rbd, dvp %$% 'par')
    kpa <- do.call(.rbd, dvp %$% 'kpa')
    
    rpt <- within(list(),
    {
        par <- DF(dat='par', mtd=rownames(par), par)
        dvp <- DF(dat='dvp', mtd=names(dvp), do.call(rbind, dvp %$% 'rpt'), rtm=rtm)
        evl <- DF(dat='evl', mtd=names(evl), do.call(.rbd, evl))
    })
    library(reshape2)
    rpt <- lapply(rpt, melt, variable.name='key', value.name='val')
    rpt <- cbind(arg, do.call(rbind, rpt))
    rownames(rpt) <- NULL

    print(list(par=par, rtm=rtm))
    rpt
}

test <- function()
{
    r <- main(N=1024, P=10000, Q=2, R=1, vcs=1, frq=.2, oks=~PL+X3)
    r <- main(N=512, P=8192, Q=2, R=1, vcs=c(0.2, 1.1, 2.0), frq=.2, oks=~LN3, seed=3)
}
