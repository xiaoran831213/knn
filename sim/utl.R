## kernels combinations
library(MASS)
source('R/utl.R')
source('sim/gsm.R')

#' get and divide a segment from the genome
#' 
#' @param gmx the genomic matrix
#' @param N the size of each cohort
#' @param P the number of variants
#' @param Q the number of training cohorts
#' @param R the number of testing cohort
get.gmx <- function(gls, N=100, P=50, Q=4, R=1)
{
    ## masks
    dvp <- c(rep(1L:Q, each=N), rep(0L, N * R))
    evl <- c(rep(0L, N * Q), rep(1L:R, each=N))
    gsp <- list()

    ## select variants and samples
    for(i in seq_along(gls))
    {
        gmx <- gls[[i]]

        ## indices
        idx <- sample.int(nrow(gmx))[seq.int(N * Q + N * R)]
        jdx <- seq(sample.int(ncol(gmx) - P, 1L), l=P)
        gmx <- gmx[idx, jdx]

        ## remove degeneracy
        af <- colMeans(gmx) / 2.0
        gmx <- gmx[, pmin(af, 1 - af) >= 0.05]
    
        ## divide
        gls[[i]] <- gmx
        gsp[[i]] <- rep(i, ncol(gmx))
    }

    ## pack up and return
    gsp <- unlist(gsp)
    gmx <- do.call(cbind, gls)
    list(gmx=gmx, dvp=dvp, evl=evl, gsp=gsp)
}

#' get variance components randomly
#'
#' @param obj the number of compoents, or a list
#' @param dfs degree of freedom for each component (rotated)
get.vcs <- function(obj, dfs=1)
{
    if(is.numeric(obj) && length(obj) == 1)
        nvc <- obj
    else
        nvc <- length(obj)
    dfs <- rep(dfs, l=nvc)
    rchisq(nvc, dfs)
}

#' get functional mask
get.fmk <- function(obj, frq=.5)
{
    if(is.numeric(obj) && length(obj) == 1)
        P <- obj
    else
        P <- NCOL(obj)
    sample(c(rep(TRUE, P * frq), rep(FALSE, P - P * frq)))
}

#' get simulated response
#'
#' @param dat genomic matrix, with training and testing masks.
#' @param frq frequency of functional variants
#' @param lnk link function to tranform the simulated signal;
#' @param eps size of noise
#' @param V the list of data generating kernels
#' @param het inter cohort heterogeneity
#' @param vc1 use these variance compoent for the shared effect,
#' the heterogeneity effect is still randomly generated.
#'
#' @return a list of lists, where the inner lists contains original genotypes and
#' generated response.
get.sim <- function(dat, frq=1, lnk=NL, eps=1, oks=~LN1, ...)
{
    dot <- list(...)
    mdl <- dot$mdl %||% A1
    svc <- dot$svc %||% 1
    vcs <- dot$vcs %||% NULL
    ep2 <- dot$ep2 %||% NULL
    vc2 <- dot$vc2 %||% NULL
    
    Q <- with(dat, length(unique(dvp[dvp > 0])))
    R <- with(dat, length(unique(evl[evl > 0])))
    one <- function(gmx, ep.=1, vc.=NULL, fmk=NULL, ...)
    {
        N <- NROW(gmx)
        within(list(gmx=gmx),
        {
            ## function kernel
            fmk <- fmk %||% get.fmk(gmx, frq)
            fmx <- gmx[, fmk]
            fmx <- gsm(mdl, fmx)
            fnl <- krn(fmx, oks)
            vc. <- vc. %||% get.vcs(fnl, svc)
            vc. <- rep(vc., l=length(fnl))
            names(vc.) <- names(fnl)
            cmx <- cmb(fnl, vc., TRUE)
            eta <- drop(mvrnorm(1, rep(0, N), cmx))
            rsp <- lnk(eta) + rnorm(N, 0, sqrt(ep.))
        })
    }
    ## core effect
    dat <- within(dat, dvp <- one(gmx[dvp > 0, ], eps, vcs))
    if(is.null(vc2))
        vc2 <- dat$dvp$vcs
    if(is.null(ep2))
        ep2 <- eps
    fmk <- dat$dvp$fmk
    dat <- within(dat, evl <- one(gmx[evl > 0, ], ep2, vc2, fmk))
    dat
}

## randomly pick a RDS file from a directory
get.rds <- function(..., n=1)
{
    ds <- list(...)
    fs <- sapply(ds, dir, "[.]rds", TRUE, TRUE)
    sample(fs, n)
}

pars <- function(dvp, ref=NULL)
{
    par <- dvp %$% 'par'
    key <- unlist(sapply(par, names), use.names=FALSE)
    if(!is.null(ref))
        key <- c(names(ref), key)
    key <- unique(key)
    mtd <- names(dvp)

    par <- sapply(par, `[`, key)
    rownames(par) <- key
    par[is.na(par)] <- 0
    par <- data.frame(t(par))
    
    if(!is.null(ref))
    {
        ref <- ref[key]
        names(ref) <- key
        ref[is.na(ref)] <- 0
        par <- rbind(par, REF=ref)
    }
    par
}

bias <- function(dvp, eps, vcs)
{
    ## bias assesment
    par <- dvp %$% 'par'                # estimates
    ref <- c(eps=eps, vcs)              # reference

    ## parameter names, and methods
    key <- unique(unlist(c(names(ref), sapply(par, names)), use.names=FALSE))
    mtd <- names(dvp)
    
    ## alignment:
    ## 1) the estimates
    par <- sapply(par, `[`, key)
    rownames(par) <- key
    par[is.na(par)] <- 0

    ## 2) reference
    ref <- ref[key]
    names(ref) <- key
    ref[is.na(ref)] <- 0

    ## the table of bias
    bia <- data.frame(t(par - ref))
    bia <- reshape(bia, key, 'val', timevar='key', idvar='mtd', ids=mtd, times=key, direction='l')
    DF(dat='bia', bia)
}
