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
source('R/agg.R')
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
    dat <- get.sim(dat, frq, lnk, eps, oks, c(rep(het, Q), rep(.0, R)))
    dvp <- dat[+(1:Q)]
    evl <- dat[-(1:Q)]
    
    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()

    ## fit each cohort separately
    sep <- list()
    dvp <- lapply(dvp, function(.) within(., knl <- krn(gmx, yks)))
    ## ******** TODO: change "fit" to director concatenation ********
    sep$mnq <- lapply(dvp, function(.) {r <- with(., knl.mnq(rsp, knl, ...)); c(., r)})
    sep$rop <- lapply(dvp, function(.) {r <- with(., rop.lmm(rsp, knl, ...)); c(., r)})
    ## sep$gct <- lapply(dvp, function(.) {r <- with(., gcta.reml(rsp, knl)); c(., r)})
    
    dvp <- with(sep,
    {
        ## combine all training samples
        gmx <- do.call(rbind, EL2(dvp, 'gmx')) # combined G
        rsp <- unlist(EL2(dvp, 'rsp'))         # combined Y
        knl <- krn(gmx, yks)                   # combined K
        rpt <- list()
        par <- list()

        ## MINQUE
        . <- mean(mnq %$% 'rpt' %[% 'val')
        rpt <- CL(rpt, DF(mtd='mnq.avg', key=rownames(.), .))
        . <- knl.mnq(rsp, knl, ...)
        rpt <- CL(rpt, DF(mtd='mnq.whl', .$rpt))
        par$mnq <- DF(mat(mnq), whl=.$par)

        ## MLE
        . <- mean(rop %$% 'rpt' %[% 'val')
        rpt <- CL(rpt, DF(mtd='mle.avg', key=rownames(.), .))
        . <- rop.lmm(rsp, knl, ...)
        rpt <- CL(rpt, DF(mtd='mle.whl', .$rpt))
        par$rop <- DF(mat(rop), whl=.$par)

        ## GCTA-REML
        ## . <- mean(gct %$% 'rpt' %[% 'val')
        ## rpt <- CL(rpt, DF(mtd='gct.avg', key=rownames(.), .))
        ## . <- gcta.reml(rsp, knl)
        ## rpt <- CL(rpt, DF(mtd='gct.whl', .$rpt))
        ## par$gct <- DF(mat(gct), whl=.$par)

        rpt <- CL(rpt, DF(mtd='nul', nul(rsp)))
        rpt <- cbind(dat='dvp', do.call(rbind, rpt))
        list(rpt=rpt, par=par)
    })
    print(dvp$par)
    ## ----------------------- Testing Errors ----------------------- ##

    evl <- with(sep,
    {
        gmx <- do.call(rbind, EL2(evl, 'gmx')) # combined G
        rsp <- unlist(EL2(evl, 'rsp'))         # combined Y
        knl <- krn(gmx, yks, ...)
        rpt <- list()

        ## MINQUE
        . <- lapply(dvp$par$mnq, vpd, y=rsp, K=knl, ...)
        . <- mapply(function(n, r) DF(mtd=paste0('mnq.', n), r), names(.), ., SIMPLIFY=FALSE)
        rpt <- c(rpt, .)
        . <- mean(lapply(mnq %$% 'par', vpd, y=rsp, K=knl, rt=0))
        rpt <- CL(rpt, DF(mtd='mnq.avg', key=names(.), val=.))

        ## MLE
        . <- lapply(dvp$par$rop, vpd, y=rsp, K=knl, ...)
        . <- mapply(function(n, r) DF(mtd=paste0('mle.', n), r), names(.), ., SIMPLIFY=FALSE)
        rpt <- c(rpt, .)
        . <- mean(lapply(rop %$% 'par', vpd, y=rsp, K=knl, rt=0))
        rpt <- CL(rpt, DF(mtd='mle.avg', key=names(.), val=.))

        ## GCTA-REML
        ## . <- lapply(dvp$par$gct, vpd, y=rsp, K=knl, ...)
        ## . <- mapply(function(n, r) DF(mtd=paste0('gct.', n), r), names(.), ., SIMPLIFY=FALSE)
        ## rpt <- c(rpt, .)
        ## . <- mean(lapply(gct %$% 'par', vpd, y=rsp, K=knl, rt=0))
        ## rpt <- CL(rpt, DF(mtd='gct.avg', key=names(.), val=.))
        
        ## NULL & GOLD
        rpt <- CL(rpt, DF(mtd='nul', nul(rsp)))
        
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
    r <- main(N=250, P=2000, Q=4, R=2, frq=.1, eps=.1, oks=ga, yks=p2, het=0.2)
}
