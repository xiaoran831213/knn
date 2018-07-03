library(magrittr)
library(Matrix)
source('R/hlp.R')
source('R/kpl.R')
source('R/gct/grm.R')
source('R/gct/gct.R')
source('R/utl.R')
source('R/lmm.R')
source("R/mnq.R")
source("R/kmq.R")
source("sim/sim_kpl.R")
source('sim/utl.R')

library(devtools)
devtools::load_all()

#' simulation of kernel deep neural network;
#' @param N size of devolpment dataset
#' @param P number of variants
#' @param H size of evaluating dataset 
#' @param frq fraction of functional variants
#' @param eps size of white noise
#' @param bsz size of mini-batches
main <- function(N, P, H=N, frq=.1, lnk=I, eps=.1, oks=c(p1), yks=c(p1), ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    arg <- match.call() %>% tail(-1) %>% as.list # %>% as.data.frame
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    gmx <- readRDS('data/p35_c05.rds')$gmx
    gmx <- sample.gmx(gmx, N, P, Q=1, H=H) 
    vcs <- c(eps=eps, vc=rchisq(length(oks), 1))
    cvs <- c(eps=id, cv=oks)
    sim <- get.sim(gmx, vcs, frq, lnk, cvs, ejt=0) # generate response
    dvp <- within(sim[[1]],
    {
        knl <- krn(gmx, yks)
        nwt <- nwt.lmm(rsp, knl)
        mnq <- kpc.mnq(rsp, knl, ...)
        rop <- rop.lmm(rsp, knl)
    })
    evl <- within(sim[[2]],
    {
        knl <- krn(gmx, yks)
        nwt <- DF(mtd='nwt', knl.prd(rsp, knl, dvp$nwt$par))
        mnq <- DF(mtd='mnq', knl.prd(rsp, knl, dvp$mnq$par, logged=FALSE, ...))
        rop <- DF(mtd='rop', knl.prd(rsp, knl, dvp$rop$par))
    })
    
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    rpt <- cl(rpt, DF(dat='dvp', mtd='mnq', dvp$mnq$rpt))
    rpt <- cl(rpt, DF(dat='dvp', mtd='rop', dvp$rop$rpt))
    rpt <- cl(rpt, DF(dat='dvp', mtd='nwt', dvp$nwt$rpt))
    rpt <- cl(rpt, DF(dat='evl', evl$mnq))
    rpt <- cl(rpt, DF(dat='evl', evl$rop))
    rpt <- cl(rpt, DF(dat='evl', evl$rpt))
    rpt <- cl(rpt, DF(dat='evl', evl$nwt))

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}

test <- function()
{
    r <- main(N=500, P=2000, H=1500, frq=.2, eps=.1, oks=c(p1), yks=c(p1))
}
