library(magrittr)
library(Matrix)
library(devtools)

source('R/hlp.R')
source('R/kpl.R')
source('R/gct/grm.R')
source('R/gct/gct.R')
source('R/utl.R')
source('R/lmm.R')
source("R/mnq.R")
source("R/kmq.R")
source("R/msg.R")
source("R/bat.R")
source("sim/sim_kpl.R")
source('sim/utl.R')

devtools::load_all()

#' simulation of kernel deep neural network;
#' @param N size of devolpment dataset
#' @param P number of variants
#' @param frq fraction of functional variants
#' @param eps size of white noise
#' @param bsz size of mini-batches
main <- function(N, P, Q=1, R=1, frq=.1, lnk=I, eps=.1, oks=p1, yks=p1, ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    het <- dot$het %||% .0
    arg <- match.call() %>% tail(-1) %>% as.list
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    gmx <- get.gmx(readRDS('data/p35_c05.rds')$gmx, N, P, Q, R)
    dat <- with(gmx, c(dvp, evl))
    dat <- get2(dat, frq, lnk, eps, oks, c(rep(het, Q), rep(het, R)))
    dvp <- dat[+(1:Q)]
    evl <- dat[-(1:Q)]
    dvp <- within(list(),
    {
        gmx <- do.call(rbind, EL2(dvp, 'gmx'))
        rsp <- unlist(EL2(dvp, 'rsp'))
        knl <- krn(gmx, yks)
        ## gct <- gcta.reml(rsp, knl)
        mnq <- knl.mnq(rsp, knl)
        kmq <- kpc.mnq(rsp, knl, ...)
        kml <- GBT(rop.lmm, rsp, knl, ...)
        kml.par <- kml$par
        kml <- do.call(rbind, lapply(kml.par, vpd, y=rsp, K=knl))
        kml.rpt <- DF(mtd=paste0('mle.', sub("[.][^.]*$", "", rownames(kml))), kml)
        mle <- rop.lmm(rsp, knl)
    })
    evl <- within(list(),
    {
        gmx <- do.call(rbind, EL2(evl, 'gmx'))
        rsp <- unlist(EL2(evl, 'rsp'))
        knl <- krn(gmx, yks)
        ## gct <- DF(mtd='gct', vpd(rsp, knl, dvp$gct$par))
        mnq <- DF(mtd='mnq.whl', vpd(rsp, knl, dvp$mnq$par))
        kmq <- DF(mtd='mnq.ssz', vpd(rsp, knl, dvp$kmq$par))
        kml <- do.call(rbind, lapply(dvp$kml.par, vpd, y=rsp, K=knl))
        kml <- DF(mtd=paste0('mle.', sub("[.][^.]*$", "", rownames(kml))), kml)
        mle <- DF(mtd='mle.whl', vpd(rsp, knl, dvp$mle$par))
    })
    
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    rpt <- CL(rpt, DF(dat='dvp', mtd='kmq', dvp$kmq$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='mnq', dvp$mnq$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='mle.whl', dvp$mle$rpt))
    rpt <- CL(rpt, DF(dat='dvp', dvp$kml.rpt))
    ## rpt <- CL(rpt, DF(dat='dvp', mtd='gct', dvp$gct$rpt))
    rpt <- CL(rpt, DF(dat='dvp', mtd='nul', nul(dvp$rsp)))

    rpt <- CL(rpt, DF(dat='evl', evl$kmq))
    rpt <- CL(rpt, DF(dat='evl', evl$mnq))
    rpt <- CL(rpt, DF(dat='evl', evl$mle))
    rpt <- CL(rpt, DF(dat='evl', evl$kml))
    ## rpt <- CL(rpt, DF(dat='evl', evl$gct))
    rpt <- CL(rpt, DF(dat='evl', mtd='nul', nul(evl$rsp)))
    
    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    print(list(kmq=dvp$kmq$par, mnq=dvp$mnq$par, gct=dvp$gct$par, mle=dvp$mle$par))
    invisible(ret)
}

test <- function()
{
    r <- main(N=250, P=3000, Q=8, R=2, frq=.1, eps=.1, oks=c(ga, p2), yks=p1, bsz=20, het=.3)
}


test.lmm <- function(N=100, P=10, t=20)
{
    library(microbenchmark)
    x <- matrix(rnorm(N * P), N, P)
    K <- krn(x, c(id, p2, ga))
    L <- length(K)
    Q <- 2
    W <- matrix(rchisq(L * Q, 1), L, Q)
    Y <- matrix(rnorm(N * Q), N, Q)

    m <- microbenchmark(
        r1 <- lmm.dv1(W, K, Y),
        r2 <- cpp.dv1(W, K, Y),
        times=t)
    print(m)
    print(all.equal(r1, r2))
    
    list(r1=r1, r2=r2)
}
