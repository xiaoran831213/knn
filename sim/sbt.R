source('R/hlp.R')                       # helpers
source('R/kpl.R')                       # kernel players
source('R/utl.R')                       # utilities
source('R/vcm.R')                       # variance component models (VCM)
source('R/msg.R')                       # message board
source('R/bat.R')                       # batched VCM trainer
source('R/agg.R')                       # model aggregation
source('R/bat.R')                       # batched trainer
source('sim/grm.R')                     # genomic relatedness matrix (plink, GCTA)
source('sim/gct.R')                     # GCTA wrapper
source('sim/utl.R')                     # simulation utilites
source('sim/mdl.R')                     # models
source('sim/lnk.R')                     # link functions
source('sim/mtd.R')                     # methods
source('sim/eps.R')                     # noise
source('sim/gen.R')

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
main <- function(N, P, Q=1, R=1, frq=.1, lnk=NL, vcs=1, oks=~L1, ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    het <- dot$het %||% 0.0
    efn <- dot$efn %||% EGS             # epsilon function for noise
    arg <- as.list(match.call()[-1])
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)

    ## ------------------------- data genration ------------------------- ##
    ## for each of the Q groups, choose N samples and P features -> Training
    dat <- sim('sim/dat', N=N, P=P, Q=Q, R=R, frq=frq, lnk=lnk, oks=oks, vcs=vcs, ...)
    dat$evl$par <- NULL
    ref <- dat$dvp$par

    ## training
    dvp <- with(dat$dvp,
    {
        ret <- list()
        kn3 <- krn(gmx, ~L3)
        ret <- CL(ret, NUL=MNQ(rsp, kn3[0:0]))
        ret <- CL(ret, GC1=GCT(rsp, kn3[1:1]))
        ret <- CL(ret, MN2=BM2(rsp, kn3[1:2]))
        ret <- CL(ret, MN3=BM2(rsp, kn3[1:3]))
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
        ret <- CL(ret, GCT=vpd(rsp, kn1, dvp$GCT$par))
        ret <- CL(ret, MNQ=vpd(rsp, kn2, dvp$MNQ$par))
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
    r <- main(N=512, P=10000, Q=2, R=1, efn=EGS, eps=2, vcs=c(2, 2), frq=.2, pss=0, bmq=ONM)
}
