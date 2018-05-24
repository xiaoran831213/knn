library(MASS)
library(magrittr)
library(Matrix)
source('R/hlp.R')
source('R/kpl.R')
source('R/gct/grm.R')
source('R/gct/gct.R')
source('R/utl.R')
source('R/lmm.R')
source("R/mnq.R")
source("sim/sim_kpl.R")

#' simulation of kernel deep neural network;
#' @param gno list contains gnomic map, subject id, and genomic matrix
#' in alternative allele dosage format;
#' @param N numeric draw this many samples for model developing;
#' @param P numeric draw this many features (i.e., SNPs);
#' @param H numeric draw this many samples for model eveluation;
#' @param frq numeric percentage of functional features (i.e., casual SNPs);
#' @param ycv character the covariance structure to generate noise free outcome;
#' l: linear kernel (cross-product); G: Gaussian; L: Laplacian
#' @param eps numeric size of white noise adds to the noise free outcome;
#' @param M numeric the number of latent features;
#' @param ukn character the basic kernel types, for the kernel network;
#' l: linear kernel (cross-product)
#' @param ikn character the inner kernel types, for the kernel network;
#' features, required to calcuate the loss and gradient;
#' @param bsz numeric batch size for mini-batch based training, when NULL or
#' greater than the sample size N, the whole data training is resumed.
main <- function(gno, N, P, H=N, frq=1, lnk=I, eps=.1, oks=c(id, p1), yks=c(id, p1), ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)
    arg <- match.call() %>% tail(-2) %>% as.list # %>% as.data.frame
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)
    
    if(is.character(gno)) gno <- readRDS(gno)
    if(is.null(gno)) gno <- readRDS('data/p35_cmn.rds')

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    gmx <- get.gmx(gno$gmx, N, H, P)
    gmx.dvp <- gmx$dvp
    gmx.evl <- gmx$evl

    ## pick functional variants
    fmk <- sample(c(rep(1, P * frq), rep(0, P - P * frq)))
    fmx.dvp <- sweep(gmx.dvp, 2L, fmk, `*`)
    fmx.evl <- sweep(gmx.evl, 2L, fmk, `*`)

    ## simulation generated variables
    sim <- within(list(),
    {
        txy <- log(rbind(phi=eps, matrix(runif(length(oks) - 1, max=3))))
        PF("PHI.SIM = %9.3f\n", exp(txy[1, ]))

        dvp <- within(list(),
        {
            ycv.fmx <- cmb(krn(fmx.dvp, oks), exp(txy))[[1]]
            ycv.gmx <- cmb(krn(gmx.dvp, oks), exp(txy))[[1]]
            y <- mvn(1, ycv.fmx) %>% drop %>% lnk
            PF("NLK.DVP.FMX = %9.3f\n", nlk(y, ycv.fmx))
            PF("NLK.DVP.GMX = %9.3f\n", nlk(y, ycv.gmx))
        })

        evl <- within(list(),
        {
            ycv.fmx <- cmb(krn(fmx.evl, oks), exp(txy))[[1]]
            ycv.gmx <- cmb(krn(gmx.evl, oks), exp(txy))[[1]]
            y <- mvn(1, ycv.fmx) %>% drop %>% lnk
            PF("NLK.EVL.FMX = %9.3f\n", nlk(y, ycv.fmx))
            PF("NLK.EVL.GMX = %9.3f\n", nlk(y, ycv.gmx))
        })
    })

    ## ----------------------- KDN Model Fitting ----------------------- ##
    rpt <- list()
    ykn.dvp <- krn(gmx.dvp, yks, ...)   # kernel for training
    ykn.evl <- krn(gmx.evl, yks, ...)   # kernel for testing
    
    ## initialize parameters
    ini.txy <- matrix(rnorm(length(yks), sd=.05), length(yks), 1)

    ## use R's optimizer
    rop.dvp <- rop.lmm(sim$dvp$y, ykn.dvp, ini.txy)
    rop.evl.rpt <- knl.prd(sim$evl$y, ykn.evl, rop.dvp$par)
    rpt <- cl(rpt, DF(mtd='rop', dat='dvp', rop.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='rop', dat='evl', rop.evl.rpt))

    ## use MINQUE
    mnq.dvp <- mnq.lmm(sim$dvp$y, ykn.dvp)
    another <- knl.mnq(sim$dvp$y, ykn.dvp)
    print(all.equal(mnq.dvp$rpt[2:4, ], another$rpt[2:4, ]))
    mnq.evl.rpt <- knl.prd(sim$evl$y, ykn.evl, mnq.dvp$par, logged=FALSE)
    rpt <- cl(rpt, DF(mtd='mnq', dat='dvp', mnq.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='mnq', dat='evl', mnq.evl.rpt))
    
    ## use GCTA:
    gct.dvp <- gcta.reml(sim$dvp$y, ykn.dvp)
    gct.evl.rpt <- knl.prd(sim$evl$y, ykn.evl, gct.dvp$par, logged=FALSE)
    rpt <- cl(rpt, DF(mtd='gct', dat='dvp', gct.dvp$rpt))
    rpt <- cl(rpt, DF(mtd='gct', dat='evl', gct.evl.rpt))

    ## oracle fit and null fit:
    rpt <- cl(rpt, DF(mtd='fmx', dat='dvp', lmm(sim$dvp$y, sim$dvp$ycv.fmx, eps)))
    rpt <- cl(rpt, DF(mtd='fmx', dat='evl', lmm(sim$evl$y, sim$evl$ycv.fmx, eps)))
    rpt <- cl(rpt, DF(mtd='gmx', dat='dvp', lmm(sim$dvp$y, sim$dvp$ycv.gmx, eps)))
    rpt <- cl(rpt, DF(mtd='gmx', dat='evl', lmm(sim$evl$y, sim$evl$ycv.gmx, eps)))
    rpt <- cl(rpt, DF(mtd='nul', dat='dvp', nul(sim$dvp$y)))
    rpt <- cl(rpt, DF(mtd='nul', dat='evl', nul(sim$evl$y)))

    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    rpt <- within(rpt, val <- round(val, 3L))
    ## ret <- list(arg=arg, rpt=rpt, sim=sim)
    ret <- cbind(arg, rpt)

    invisible(ret)
}
