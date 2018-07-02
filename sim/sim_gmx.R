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
#' @param gno gnomic map, subject id, and genomic matrix in dosage format;
#' @param N numeric draw this many samples for model developing;
#' @param P numeric draw this many features (i.e., SNPs);
#' @param H numeric draw this many samples for model eveluation;
#' @param frq numeric percentage of functional features (i.e., casual SNPs);
#' @param eps numeric size of white noise adds to the noise free outcome;
#' @param bsz numeric batch size for mini-batch based training, when NULL or
main <- function(gno, N, P, H=N, frq=.1, lnk=I, eps=.2, oks=c(id, p1), yks=c(id, p1), ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    arg <- match.call() %>% tail(-2) %>% as.list # %>% as.data.frame
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)

    if(is.character(gno)) gno <- readRDS(gno)
    if(is.null(gno)) gno <- readRDS('data/p35_c05.rds')

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    nvc <- length(oks)
    gmx <- sample.gmx(gno$gmx, N, P, Q=2, H=H) 
    vcs <- c(eps, rchisq(nvc - 1, 1))

    sim <- get.sim(gmx, vcs, frq, lnk, oks, ejt=0)
    dvp <- within(sim[[1]],
    {
        knl <- krn(gmx, yks)
        mnq <- kpc.mnq(rsp, knl[-1], ...)
        rop <- rop.lmm(rsp, knl)
        nwt <- nwt.lmm(rsp, knl)
    })
    evl <- within(sim[[2]],
    {
        knl <- krn(gmx, yks)
        mnq <- DF(mtd='mnq', knl.mnq.evl(rsp, knl[-1], dvp$mnq$par))
        rop <- DF(mtd='rop', knl.prd(rsp, knl, dvp$rop$par))
        nwt <- DF(mtd='nwt', knl.prd(rsp, knl, dvp$nwt$par))
    })
    
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
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
    r <- main(NULL, N=500, P=3000, H=1500, frq=.1,  eps=.1, oks=c(id, p1), yks=c(id, p1),
              ejt=0.0, bsz=100, wep=2)
}
