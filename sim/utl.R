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
get.sim <- function(dat, frq=1, lnk=i1, eps=1, oks=~p1, yks=oks, ...)
{
    dot <- list(...)
    ## mdl <- dot$mdl %||% a1
    svc <- dot$svc %||% 1
    
    Q <- with(dat, length(unique(dvp[dvp > 0])))
    R <- with(dat, length(unique(evl[evl > 0])))
    one <- function(gmx, vcs=NULL, fmk=NULL, ...)
    {
        P <- with(dat, NCOL(gmx))
        N <- with(dat, NROW(gmx))

        within(list(gmx=gmx),
        {
            ## working kernel
            knl <- krn(gmx, yks)        

            ## function kernel
            fmk <- fmk %||% get.fmk(gmx, frq)
            ## fmx <- gsm(mdl, gmx[, fmk], rm.dup=FALSE)
            fmx <- gmx[, fmk]
            fnl <- krn(fmx, oks)
            vcs <- vcs %||% get.vcs(fnl, svc)
            cmx <- cmb(fnl, vcs)[[1]] #  + diag(rnorm(nrow(gmx), 0, 1e-4))
            eta <- drop(mvrnorm(1, rep(0, nrow(gmx)), cmx))
            nos <- rnorm(nrow(gmx), 0, sqrt(eps))
            rsp <- lnk(eta + nos)
        })
    }
    ## core effect
    dat <- within(dat, dvp <- one(gmx[dvp > 0, ]))
    vcs <- dat$dvp$vcs
    fmk <- dat$dvp$fmk
    dat <- within(dat, evl <- one(gmx[evl > 0, ], vcs, fmk))
    dat
}

## randomly pick a RDS file from a directory
get.rds <- function(..., n=1)
{
    ds <- list(...)
    fs <- sapply(ds, dir, "[.]rds", TRUE, TRUE)
    sample(fs, n)
}

bias <- function(dvp, eps, vcs)
{
    ## bias assesment
    par <- dvp %$% 'par'                # estimates
    ref <- c(eps, vcs)                  # reference
    len <- max(length(ref), unlist(par %$% length))
    key <- paste0('bia.', c('eps', paste0('vc', 1:(len-1))))

    ## zero padding
    par <- lapply(par, function(.) c(., rep(0, len - length(.))))
    ref <- c(ref, rep(0, len - length(ref)))

    ## bias
    ret <- lapply(names(par), function(.)
    {
        DF(mtd=., key=key, val=par[[.]] - ref)
    })
    DF(dat='dvp', do.call(rbind, ret))
}
