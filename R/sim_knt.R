library(dplyr, warn.conflicts=FALSE)
library(rstiefel)
source('R/hlp.R')
source('R/kpl1.R')
source('R/kd3.R')
#' simulation of kernel deep neural network;
#' @param gno list contains gnomic map, subject id, and genomic matrix
#' in alternative allele dosage format;
#' @param N numeric draw this many samples for model developing;
#' @param P numeric draw this many features (i.e., SNPs);
#' @param H numeric draw this many samples for model eveluation;
#' @param frq numeric percentage of functional features (i.e., casual SNPs);
#' @param ycv character the covariance structure to generate noise free outcome;
#' l: linear kernel (cross-product); G: Gaussian; L: Laplacian
#' @param PHI numeric size of white noise adds to the noise free outcome;
#' @param M numeric the number of latent features;
#' @param ukn character the basic kernel types, for the kernel network;
#' l: linear kernel (cross-product)
#' @param ikn character the inner kernel types, for the kernel network;
#' features, required to calcuate the loss and gradient;
#' @param bsz numeric batch size for mini-batch based training, when NULL or
#' greater than the sample size N, the whole data training is resumed.
#' @param mmt numeric the mommentum of stochastical gradient descent; the actual
#' direction of parameter adjustment is the weight sum of previous direction and
#' the current gradient, and the weight of previous step is the mommentum;
#' @param max.itr integer the maximum allowed traing steps;
sm3 <- function(gno, N, P, H=N, frq=.5, lnk=I, PHI=1, M=10, ukn='p', ykn='p', ...)
{
    options(stringsAsFactors=FALSE)
    dot <- list(...)
    set.seed(dot$seed)

    .ln <- substitute(lnk) %>% as.character %>% paste(collapse='')
    arg <- match.call() %>% tail(-2) %>% as.list # %>% as.data.frame
    arg <- within(arg,
    {
        lnk <- .ln
        ukn <- knl2str(eval(ukn))             # basic kernel
        ykn <- knl2str(eval(ykn))             # inner kernel
    }) %>% as.data.frame
    
    if(is.character(gno)) gno <- readRDS(gno)
    if(is.null(gno)) gno <- readRDS('data/p35_005.rds')

    ## ------------------------- data genration ------------------------- ##
    ## choose N samples and P features for both development and evaluation
    gmx <- gno$gmx
    N <- min(N, nrow(gmx) - H)
    P <- min(P, ncol(gmx))
    idx <- sample.int(nrow(gmx), N + H) %>% sort
    jdx <- sample.int(ncol(gmx) - P, 1) %>% seq(l=P)
    gmx.dvp <- as.matrix(gmx[idx[+(1:N)], jdx])
    gmx.evl <- as.matrix(gmx[idx[-(1:N)], jdx])
    gmx.pmu <- as.matrix(gmx[sample(1:nrow(gmx), N), sample(1:ncol(gmx), P)])

    ## pick functional variants
    fmk <- sample(c(rep(1, P * frq), rep(0, P - P * frq)))
    fmx.dvp <- sweep(gmx.dvp, 2L, fmk, `*`)
    fmx.evl <- sweep(gmx.evl, 2L, fmk, `*`)
    
    ## simulation generated variables
    L <- length(ukn)
    J <- length(ykn)
    txy.sim <- log(rbind(phi=PHI, matrix(1, J - 1, 1)))
    PF("PHI.SIM = %9.3f\n", exp(txy.sim[1, ]))

    ycv.dvp <- cmb(findKernel(fmx.dvp, ykn), exp(txy.sim))[[1]]
    ych.dvp <- chol(ycv.dvp)
    yac.dvp <- chol2inv(ych.dvp)
    y.dvp <- mvn(1, ych.dvp, 1) %>% drop %>% lnk
    PF("NLK.FMX = %9.3f\n", nlk(y.dvp, ycv.dvp, ych.dvp, yac.dvp))

    ycv.gmx <- cmb(findKernel(gmx.dvp, ykn), exp(txy.sim))[[1]]
    ych.gmx <- chol(ycv.gmx)
    yac.gmx <- chol2inv(ych.gmx)
    PF("NLK.GMX = %9.3f\n", nlk(y.dvp, ycv.gmx, ych.gmx, yac.gmx))

    ycv.evl <- cmb(findKernel(fmx.evl, ykn), exp(txy.sim))[[1]]
    ych.evl <- chol(ycv.evl)
    yac.evl <- chol2inv(ych.evl)
    y.evl <- mvn(1, ych.evl, 1) %>% drop %>% lnk
    PF("NLK.EVL = %9.3f\n", nlk(y.evl, ycv.evl, ych.evl, yac.evl))

    y.pmu <- sample(c(y.dvp, y.evl), N)
    ## ----------------------- KDN model fitting ----------------------- ##
    rpt <- list()
    ## calculate basic kernels
    ukn.dvp <- findKernel(gmx.dvp, ukn, ...)
    ukn.pmu <- findKernel(gmx.pmu, ukn, ...)
    ukn.evl <- findKernel(gmx.evl, ukn, ...)

    ## initialize fitted parameters
    par.dvp <- list(tuy=matrix(rnorm(J * 1, sd=.5), J, 1),
                    txu=matrix(rnorm(L * M, sd=.5), L, M))

    ## initialize context
    ctx.dvp <- list(par=par.dvp, ukn=ukn.dvp, ykn=ykn, ycv.ref=ycv.dvp, y=y.dvp, M=M)
    ctx.pmu <- within(ctx.dvp, {ukn <- ukn.pmu; y=y.pmu})
    ctx.evl <- within(ctx.dvp, {ukn <- ukn.evl; y=y.evl; ycv.ref=ycv.evl})
    
    ## Benchmark of non KDN methods
    cfg <- select(arg, -M, -ukn, -ykn)
    rpt <- cl(rpt, df(cfg, mtd='lmm', dat='dvp', lmm(y.dvp, ycv.dvp, PHI)))
    rpt <- cl(rpt, df(cfg, mtd='lmm', dat='evl', lmm(y.evl, ycv.evl, PHI)))
    rpt <- cl(rpt, df(cfg, mtd='lmm', dat='gmx', lmm(y.dvp, ycv.gmx, PHI)))
    rpt <- cl(rpt, df(cfg, mtd='nul', dat='dvp', nul(y.dvp)))
    ## ---------------------------------------------------------------------- ##    
    
    ## KDN Gradient Descent Now!
    fit <- sg3(ctx.dvp, ...)
    pmu <- up3(ctx.pmu, -1, ...)
    evl <- up3(ctx.evl, -1, ...)
    hst <- within(fit$hst, cbind(arg, stt))

    ## collect statistics.
    print('KDN-DVP')
    cfg <- mutate(arg, mtd='kdn', dat='dvp')
    .ky <- c('lxy', 'mse', 'loo', 'cyh')
    rpt <- cl(rpt, df(cfg, key=.ky, val=unlist(fit[.ky])))

    ## ---------------------- Evaluation ------------------------- ##
    ## report and return
    rpt <- Reduce(function(a, b) merge(a, b, all=TRUE), rpt)
    par <- list(fit=fit$par, dvp=par.dvp, sim=txy.sim)
    rpt <- within(rpt, val <- round(val, 2L))
    ret <- list(rpt=rpt, hst=hst, par=par, fit=fit, lw1=fit$lw1, lw2=fit$lw2, pmu=pmu, evl=evl)
    invisible(ret)
}
