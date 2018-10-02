## kernels combinations
library(MASS)
source('R/kpl.R')
source('R/utl.R')

ID <- function(x) list(ID=diag(1, NROW(x)))
JX <- function(x) list(JX=matrix(1, NROW(x), NROW(x)))
LN <- function(x) list(LN=ply(scale(x), degree=1))
PL <- function(x) list(PL=ply(x, degree=1))
LP <- function(x) list(LP=lap(x))
GS <- function(x) list(GS=gau(x))
IS <- function(x) list(IS=ibs(x))
KN <- function(x) list(KN=kin(x))
RS <- function(x) list(RS=ply((x >  1) / 2 - 1))
DM <- function(x) list(DM=ply((x >  0) / 2 - 1))
AD <- function(x) list(AD=ply((x     ) / 2 - 1))
HT <- function(x) list(HT=ply((x == 1) / 2 - 1))

#' polynomial expansion by kernel
#'
#' ...: the series of basic kernel functions
PK <- function(..., d=2, orth=FALSE, J=FALSE)
{
    env <- environment()
    function(x)
    {
        print(env[['coef']])
        N <- NROW(x)                    # sample size
        mtx <- matrix(0, N, N)          # sample matrix
        utr <- upper.tri(mtx, TRUE)     # sample upper.tri

        bas <- c(...)                   # bases
        bas <- unlist(lapply(bas, do.call, list(x=x)), FALSE)
        bas <- lapply(bas, `[`, utr)
        nms <- names(bas)
        
        ## expansion
        arg <- c(bas, list(degree=d, coefs=env[['coef']], raw=!orth))
        bas <- do.call(polym, arg)
        
        ## keep orthgonal coeficients
        cfs <- attr(bas, 'coefs')
        env[['coef']] <- if(length(nms) >1) cfs else if(is.null(cfs)) NULL else list(cfs)
        print(env[['coef']])

        bas <- as.data.frame(bas)
        colnames(bas) <- lapply(strsplit(colnames(bas), '[.]'), function(u)
        {
            paste0(nms, u, collapse='.')
        })
        
        lapply(bas, function(k)
        {
            mtx[utr] <- k; mtx <- t(mtx); mtx[utr] <- k; mtx
        })
    }
}

OL2 <- PK(LN, d=2, orth=TRUE)
OL3 <- PK(LN, d=3, orth=TRUE)
RL2 <- PK(LN, d=2, orth=FALSE)
RL3 <- PK(LN, d=2, orth=FALSE)
    
PG1 <- PK(PL, GS, d=1, orth=FALSE)
PG2 <- PK(PL, GS, d=2, orth=FALSE)
PG3 <- PK(PL, GS, d=3, orth=FALSE)


coef.function <- function(k)
{
    environment(k)[['coef']]
}
`coef<-` <- function(x, value)
{
    env <- environment(x)
    if(!is.null(env))
        env[['coef']] <- value
    x
}

## link functions
#' Distribution Convertor
#'
#' covert an assumed empirical nromal distribution to another distribution.
#'
#' x: the data that descript the empirical normal
#' d: the target distribution
#' .: additional paramters required by the target quantile (e.g., degree
#' of freedom for t and chisq, min and min for unif, etc.)
dc <- function(x, d=NULL, curb=0.01, ...)
{
    v <- var(x)                         # variance
    s <- sd(x)                          # sd
    m <- mean(x)                        # mean
    n <- length(x)                      # numbers
    p <- pnorm(x, m, s)                 # quantile

    p <- curb / 2 + p * (1 - curb)

    ## mean of target distribution.
    y <- switch(d, bn={2 * v}, ps={v}, ex={1 / s}, ch={v}, 0)

    z <- switch(
        d,
        ## st={     qt(p, (2 * v) / max(v - 1, 0))},
        st={     qt(p, 1 / s)},
        bn={ qbinom(p, ceiling(4 * v), .5)},
        ps={  qpois(p, v)},
        ex={   qexp(p, 1 / s)},
        ## ca={qcauchy(p) -> .; s * . / sd(.)},
        ca={qcauchy(p, 0, s)},
        ch={ qchisq(p, v / 2)})

    ## centered?
    cdc <- as.logical(get0('cdc', ifnotfound=FALSE))
    print(list(cdc=cdc))
    if(cdc)
        z <- z - y
    z
}

dc.test <- function()
{
    x <- rnorm(50000, 0, 1.5)
    x.st <- dc(x, 'st', 0.10)
    x.ca <- dc(x, 'ca', 0.05)
    print(list(nm=summary(x), st=summary(x.st), ca=summary(x.ca)))

    d <- rbind(
        DF(dst='nm', val=x),
        DF(dst='st', val=x.st),
        DF(dst='ca', val=x.ca))

    library(ggplot2)
    g <- ggplot(d, aes(val))
    g <- g + geom_freqpoly(aes(color=dst), binwidth = 0.1)
    g <- g + xlim(-15, 15)
    ## g <- g + ylim(c(1, NA))
    g
}

ST <- function(x) dc(x, 'st', 0.05)
BN <- function(x) dc(x, 'bn')
PS <- function(x) dc(x, 'ps')
EX <- function(x) dc(x, 'ex')
CA <- function(x) dc(x, 'ca', 0.05)
CH <- function(x) dc(x, 'ch')
I1 <- function(x) x
I2 <- function(x) drop(scale((x + 1)^2, scale=FALSE))
I3 <- function(x) drop(scale((x + 1)^3, scale=FALSE))
SN <- function(x) sin(2 * pi * x)


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
