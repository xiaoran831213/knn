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
CT <- function(k, s=1)
{
    k <- k - outer(rowMeans(k), colMeans(k), `+`) + mean(k)
    if(s)
        k <- k / mean(diag(k))
    k
}

LN <- function(x, o=1)
{
    k <- tcrossprod(scale(x)) / NCOL(x)
    SEQ <- c(LN1=1, LN2=2, LN3=3, LN4=4, LN5=5)[seq(2, l=o-1)]
    LNX <- lapply(SEQ, function(i) CT(k^i))
    c(list(LN1=k), LNX)
}
L1 <- function(x) LN(x, 1)
L2 <- function(x) LN(x, 2)
L3 <- function(x) LN(x, 3)
I2 <- function(x) list(I2=CT(pqw(scale(x), q=2)))
I3 <- function(x) list(I3=CT(pqw(scale(x), q=3)))

GS <- function(x) list(GS1=CT(gau(x), 0))
LP <- function(x) list(LP1=CT(lap(scale(x)), 1))
XK <- function(x)
{
    c(L3(x), GS(x))
}

FWD <- function(rsp, kns, xmx=NULL, ...) fwd(rsp, kns, xmx, tol=1e-4, rpt=1, ...)
MNQ <- function(rsp, kns, xmx=NULL, ...) mnq(rsp, kns, xmx, tol=1e-4, rpt=1, ...)

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
        ret <- CL(ret, BM2=BM2(rsp, kn3[1:3], tlr=0.05/2))
        ret <- CL(ret, BM3=BM3(rsp, kn3[1:3], tlr=0.05/2))
        ret <- CL(ret, BM4=BM4(rsp, kn3[1:3], tlr=0.05/2))
        ret
    })

    ## bias assesment
    bpa <- dvp %$% 'bpa'
    par <- do.call(rbd, dvp %$% 'par')
    par <- rbd(REF=ref, par)
    rtm <- do.call(rbd, dvp %$% 'rtm')
    dvp <- do.call(rbd, dvp %$% 'rpt')
    
    ## testing
    evl <- with(dat$evl,
    {
        ret <- list()
        kn3 <- krn(gmx, ~L3)
        ret <- CL(ret, NUL=vpd(rsp, kn3[0:0], w=par["NUL", ]))
        ret <- CL(ret, GC1=vpd(rsp, kn3[1:1], w=par["GC1", ]))
        ret <- CL(ret, BM2=vpd(rsp, kn3[1:3], w=par["BM2", ]))
        ret <- CL(ret, BM3=vpd(rsp, kn3[1:3], w=par["BM3", ]))
        ret <- CL(ret, BM4=vpd(rsp, kn3[1:3], w=par["BM4", ]))
        ret
    })
    evl <- do.call(rbd, evl)

    ## ----------------------- generate reports ----------------------- ##
    rpt <- within(list(),
    {
        par <- DF(dat='par', mtd=rownames(par), par)
        rtm <- DF(dat='rtm', mtd=rownames(rtm), rtm)
        dvp <- DF(dat='dvp', mtd=rownames(dvp), dvp)
        evl <- DF(dat='evl', mtd=rownames(evl), evl)
    })
    library(reshape2)
    rpt <- lapply(rpt, melt, id.var=c('dat', 'mtd'), variable.name='key', value.name='val')
    rpt <- cbind(arg, do.call(rbind, rpt))
    rownames(rpt) <- NULL

    print(list(par=par, rtm=rtm))
    rpt
}

test <- function()
{
    r=main(N=512, P=8192, Q=2, R=2, frq=.10, lnk=NL, oks=~L2, vcs=c(1.0, 0.0, 1.0, 0.0), seed=4)
}
