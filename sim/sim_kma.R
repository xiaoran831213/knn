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
source('R/msg.R')
source("sim/sim_kpl.R")
source('sim/utl.R')

library(devtools)
devtools::load_all()

#' simulation of kernel deep neural network;
#' @param gno gnomic map, subject id, and genomic matrix in dosage format;
#' @param N draw this many samples for training and test, for each cohort;
#' @param P draw this many features (i.e., SNPs);
#' @param frq numeric percentage of functional features (i.e., casual SNPs);
#' @param eps numeric size of white noise adds to the noise free outcome;
#' @param bsz numeric batch size for mini-batch based training, when NULL or
main <- function(N, P, Q=5, R=1, frq=.1, lnk=I, eps=.2, oks=p1, yks=p1, ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    het <- dot$het %||% .2
    arg <- match.call() %>% tail(-1) %>% as.list # %>% as.data.frame
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)
    
    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    gmx <- readRDS('data/p35_c05.rds')$gmx
    gmx <- get.gmx(gmx, N, P, Q, R)
    dat <- with(gmx, c(dvp, evl))
    dat <- get2(dat, frq, lnk, eps, oks, c(rep(het, Q), rep(.0, R)))
    dvp <- dat[+(1:Q)]
    evl <- dat[-(1:Q)]
    
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()

    ## fit each cohort separately
    sep <- list()
    dvp <- lapply(dvp, function(.) within(., knl <- krn(gmx, yks)))
    ## ******** TODO: change "fit" to director concatenation ********
    sep$mnq <- lapply(dvp, function(.) within(., fit <- knl.mnq(rsp, knl, ...)))
    sep$rop <- lapply(dvp, function(.) within(., fit <- rop.lmm(rsp, knl, ...)))
    ## sep$gct <- lapply(dvp, function(.) within(., fit <- gcta.reml(rsp, knl)))
    ## sep$kmq <- lapply(dvp, function(.) within(., fit <- kpc.mnq(rsp, knl, ...)))
    
    dvp <- with(sep,
    {
        ## combine all training samples
        gmx <- do.call(rbind, EL2(dvp, 'gmx')) # combined G
        rsp <- unlist(EL2(dvp, 'rsp'))  # combined Y
        knl <- krn(gmx, yks)            # combined K
        rpt <- list()
        par <- list()

        tmp <- mean(sapply(mnq, function(.) subset(.$fit$rpt, key=='nlk')$val))
        rpt <- CL(rpt, DF(mtd='avg.mnq', key='nlk', val=tmp))
        mnq <- knl.mnq(rsp, knl, ...)
        rpt <- CL(rpt, DF(mtd='whl.mnq', mnq$rpt))
        par$mnq <- mnq$par
        
        tmp <- mean(sapply(rop, function(.) subset(.$fit$rpt, key=='nlk')$val))
        rpt <- CL(rpt, DF(mtd='avg.rop', key='nlk', val=tmp))
        rop <- rop.lmm(rsp, knl, ...)
        rpt <- CL(rpt, DF(mtd='whl.rop', rop$rpt))
        par$rop <- rop$par
        
        ## gct <- gcta.reml(rsp, knl)
        ## rpt <- CL(rpt, DF(mtd='whl.gct', gct$rpt))
        ## par$gct <- gct$par

        ## kmq <- kpc.mnq(rsp, knl, ...)
        ## rpt <- CL(rpt, DF(mtd='whl.kmq', kmq$rpt))
        ## par$kmq <- kmq$par

        rpt <- CL(rpt, DF(mtd='nul', nul(rsp)))
        rpt <- cbind(dat='dvp', do.call(rbind, rpt))
        list(rpt=rpt, par=par)
    })

    ## ----------------------- Testing Errors ----------------------- ##
    agg <- lapply(sep, mat)
    evl <- with(sep,
    {
        gmx <- do.call(rbind, EL2(evl, 'gmx')) # combined G
        rsp <- unlist(EL2(evl, 'rsp'))         # combined Y
        knl <- krn(gmx, yks, ...)
        rpt <- list()

        ## performance: meta analysis
        rpt <- CL(rpt, DF(mtd='nlk.mnq', vpd(rsp, knl, agg$mnq[, 'nlk'], ...)))
        rpt <- CL(rpt, DF(mtd='ssz.mnq', vpd(rsp, knl, agg$mnq[, 'ssz'], ...)))
        rpt <- sapply(mnq, function(.) vpd(rsp, knl, .$fit$par, rt=0)) %>%
            rowMeans %>% {DF(mtd='avg.mnq', key=names(.), val=.)} %>% {CL(rpt, .)}

        rpt <- CL(rpt, DF(mtd='nlk.rop', vpd(rsp, knl, agg$rop[, 'nlk'], ...)))
        rpt <- CL(rpt, DF(mtd='ssz.rop', vpd(rsp, knl, agg$rop[, 'ssz'], ...)))
        rpt <- sapply(rop, function(.) vpd(rsp, knl, .$fit$par, rt=0)) %>%
            rowMeans %>% {DF(mtd='avg.rop', key=names(.), val=.)} %>% {CL(rpt, .)}

        ## rpt <- CL(rpt, DF(mtd='nlk.gct', vpd(rsp, knl, agg$gct[, 'nlk'], ...)))
        ## rpt <- CL(rpt, DF(mtd='ssz.gct', vpd(rsp, knl, agg$gct[, 'ssz'], ...)))
        ## rpt <- sapply(gct, function(.) vpd(rsp, knl, .$fit$par, rt=0)) %>%
        ##     rowMeans %>% {DF(mtd='avg.gct', key=names(.), val=.)} %>% {CL(rpt, .)}

        ## rpt <- CL(rpt, DF(mtd='nlk.kmq', vpd(rsp, knl, agg$kmq[, 'nlk'], ...)))
        ## rpt <- CL(rpt, DF(mtd='ssz.kmq', vpd(rsp, knl, agg$kmq[, 'ssz'], ...)))
        ## rpt <- sapply(kmq, function(.) vpd(rsp, knl, .$fit$par, rt=0)) %>%
        ##      rowMeans %>% {DF(mtd='avg.kmq', key=names(.), val=.)} %>% {CL(rpt, .)}
        
        ## performance: model of whole training data
        rpt <- CL(rpt, DF(mtd='whl.mnq', vpd(rsp, knl, dvp$par$mnq, ...)))
        rpt <- CL(rpt, DF(mtd='whl.rop', vpd(rsp, knl, dvp$par$rop, ...)))
        ## rpt <- CL(rpt, DF(mtd='whl.gct', vpd(rsp, knl, dvp$par$gct, ...)))
        ## rpt <- CL(rpt, DF(mtd='whl.kmq', vpd(rsp, knl, dvp$par$kmq, ...)))

        ## NULL & GOLD
        rpt <- CL(rpt, DF(mtd='nul', nul(rsp)))
        ## rpt <- CL(rpt, evl[[1]]$rpt)
        
        ## report
        rpt=cbind(dat='evl', do.call(rbind, rpt))
        list(rpt=rpt)
    })
    rpt <- list(dvp=dvp$rpt, evl=evl$rpt)

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 4L))
    ret <- cbind(arg, rpt)

    invisible(ret)
}


test <- function()
{
    r <- main(N=200, P=2000, Q=5, R=1, frq=.1, eps=.1, oks=ga, yks=p2, het=0.2)
}
