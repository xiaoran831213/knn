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

library(devtools)                       # enable the C++ functions
devtools::load_all()

#' simulation of kernel deep neural network;
#' @param N size of a sub-population
#' @param P number of variants
#' @param Q number of developing population
#' @param R number of evaluation population
#' @param frq fraction of functional variants
#' @param lnk link function to transform the raw response
#' @param eps size of white noise
#' @param oks kernels for simulated data
#' @param yks kernels for the model to be tested.
#' @param bsz size of mini-batches
#' @param sc effect scaler
main <- function(N, P, Q=1, R=1, frq=.1, lnk=I, eps=.1, oks=p1, yks=p1, ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    het <- dot$het %||% .0
    sc <- dot$sc %||% 1
    arg <- match.call() %>% tail(-1) %>% as.list
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)
    gmx <- dot$gds %||% 'data/1kg_c05.rds'

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for Q development
    ## cohorts and and R evaluation cohorts
    gmx <- get.gmx(readRDS(gmx), N, P, Q, R)
    dat <- with(gmx, c(dvp, evl))
    dat <- get.sim(dat, frq, lnk, eps, oks, het=c(rep(het, Q), rep(het, R)), vc1=NULL, sc=sc)
    dvp <- dat[+(1:Q)]
    evl <- dat[-(1:Q)]
    dvp <- within(list(),
    {
        gmx <- do.call(rbind, EL2(dvp, 'gmx')) # stacked genomics
        rsp <- unlist(EL2(dvp, 'rsp'))         # stacked response
        knl <- krn(gmx, yks)                   # kernels from the entire data
        gct <- gcta.reml(rsp, krn(gmx, p1))    # whole sample GCTA-REML
        mnq <- knl.mnq(rsp, knl, cpp=TRUE)     # whole sample MINQUE
        mle <- rop.vcm(rsp, knl, cpp=TRUE)     # whold sample MLE

        ## batched MINQUE
        kmq <- GBT(knl.mnq, rsp, knl, cpp=TRUE, ...)
        kmq.par <- kmq$par
        kmq.rtm <- kmq$rtm
        kmq <- do.call(rbind, lapply(kmq.par, vpd, y=rsp, K=knl))
        kmq.rpt <- DF(mtd=paste0('mnq.', sub("[.][^.]*$", "", rownames(kmq))), kmq)
        kmq.rpt <- rbind(kmq.rpt, DF(mtd='mnq.bat', key='rtm', val=kmq.rtm))

        ## batched MLE
        kml <- GBT(rop.vcm, rsp, knl, cpp=TRUE, ...)
        kml.par <- kml$par
        kml.rtm <- kml$rtm
        kml <- do.call(rbind, lapply(kml.par, vpd, y=rsp, K=knl))
        kml.rpt <- DF(mtd=paste0('mle.', sub("[.][^.]*$", "", rownames(kml))), kml)
        kml.rpt <- rbind(kml.rpt, DF(mtd='mle.bat', key='rtm', val=kml.rtm))
    })
    evl <- within(list(),
    {
        gmx <- do.call(rbind, EL2(evl, 'gmx')) # stacked genomics
        rsp <- unlist(EL2(evl, 'rsp'))         # stacked response
        knl <- krn(gmx, yks)                   # kernels from the entire data

        gct <- DF(mtd='gct', vpd(rsp, krn(gmx, p1), dvp$gct$par)) # evaluate GCTA-REML
        mnq <- DF(mtd='mnq.whl', vpd(rsp, knl, dvp$mnq$par))      # evaluate MINQUE
        mle <- DF(mtd='mle.whl', vpd(rsp, knl, dvp$mle$par))      # evaluate MLE

        kmq <- do.call(rbind, lapply(dvp$kmq.par, vpd, y=rsp, K=knl))
        kmq <- DF(mtd=paste0('mnq.', sub("[.][^.]*$", "", rownames(kmq))), kmq)
        
        kml <- do.call(rbind, lapply(dvp$kml.par, vpd, y=rsp, K=knl))
        kml <- DF(mtd=paste0('mle.', sub("[.][^.]*$", "", rownames(kml))), kml)
    })
    
    ## ----------------------- generate reports ----------------------- ##
    rpt <- list()
    rpt <- CL(rpt, DF(dat='dvp', mtd='mnq.whl', dvp$mnq$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='mle.whl', dvp$mle$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='gct', dvp$gct$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='nul', nul(dvp$rsp)))

    rpt <- CL(rpt, DF(dat='evl', evl$mnq))
    rpt <- CL(rpt, DF(dat='evl', evl$mle))
    rpt <- CL(rpt, DF(dat='evl', evl$gct))
    rpt <- CL(rpt, DF(dat='evl', mtd='nul', nul(evl$rsp)))

    rpt <- CL(rpt, DF(dat='dvp', dvp$kmq.rpt))
    rpt <- CL(rpt, DF(dat='evl', evl$kmq))

    rpt <- CL(rpt, DF(dat='dvp', dvp$kml.rpt))
    rpt <- CL(rpt, DF(dat='evl', evl$kml))
    
    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    print(list(kmq=dvp$kmq$par, mnq=dvp$mnq$par, gct=dvp$gct$par, mle=dvp$mle$par))
    ## print(list(mnq=dvp$mnq$par, gct=dvp$gct$par, mle=dvp$mle$par))
    invisible(ret)
}

test <- function()
{
    r <- main(N=100, P=4000, Q=8, R=15, frq=.01, eps=.1, lnk=i2, oks=p1, yks=c(ga, p2), bsz=200, sc=5)
}


test.vcm <- function(N=100, P=10, t=20)
{
    library(microbenchmark)
    x <- matrix(rnorm(N * P), N, P)
    K <- krn(x, c(id, p2, ga))
    L <- length(K)
    Q <- 2
    W <- matrix(rchisq(L * Q, 1), L, Q)
    Y <- matrix(rnorm(N * Q), N, Q)

    m <- microbenchmark(
        r1 <- vcm.dv1(W, K, Y),
        r2 <- cpp.dv1(W, K, Y),
        times=t)
    print(m)
    print(all.equal(r1, r2))
    
    list(r1=r1, r2=r2)
}
