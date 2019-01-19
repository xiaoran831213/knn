## kernels combinations
library(MASS)
source('R/utl.R')

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

#' get simulated response for once
#'
#' @param gmx N x P genomic data matrix
#' @param oks formula of output kernels
#' @param eps noise size, or functional
#' @param vcs 
get.one <- function(gmx, oks=~LN, lnk=NL, frq=.1, eps=1, ...)
{
    N <- nrow(gmx)
    dot <- list(...)

    ## functional variants
    fmk <- dot$fmk %||% get.fmk(gmx, frq) # functional mask
    fnl <- krn(gmx[, fmk], oks)           # functional kernel

    ## Variance components
    vcs <- dot$vcs %||% 1                 #
    vcs <- rep(vcs, l=length(fnl))        # propagate
    names(vcs) <- names(fnl)

    ## covariance matrix
    cmx <- cmb(fnl, vcs, TRUE)

    ## noise function
    efn <- dot$efn %||% EGS

    ## transform and return
    eta <- drop(mvrnorm(1, rep(0, N), cmx))
    rsp <- lnk(eta) + efn(N, eps)

    ret <- list(gmx=gmx, rsp=rsp, fmk=fmk, vcs=vcs)
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
get.sim <- function(dat, ...)
{
    dot <- list(...)
    Q <- with(dat, length(unique(dvp[dvp > 0]))) # dvp groups
    R <- with(dat, length(unique(evl[evl > 0]))) # evl groups

    ## core effect
    dat <- within(dat, dvp <- get.one(gmx[dvp > 0, ], ...))
    dat <- within(dat, evl <- get.one(gmx[evl > 0, ], ..., fmk=dvp$fmk))
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
    key <- unlist(lapply(par, names), use.names=FALSE)
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
    ref <- c(EPS=eps, vcs)              # reference

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
