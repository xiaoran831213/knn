## kernels combinations
library(MASS)
source('R/kpl.R')
id <- c(id=function(x) diag(1, NROW(x)))
p1 <- c(o1=function(x) ply(x, degree=1))
p2 <- c(o2=function(x) ply(x, degree=2))
ga <- c(ga=function(x) gau(x))
lp <- c(lp=function(x) lap(x))
ib <- c(ib=function(x) ibs(x))

gk <- c(ga=function(x) gau(x),
        ki=function(x) kin(x))
kg <- c(ki=function(x) kin(x),
        ga=function(x) gau(x))

k1 <- c(k1=function(x) kin(x))          # kinship
k2 <- c(k2=function(x) kin(x)^2)        # kinship^2
s1 <- c(s1=function(x) esn(x, p=1.0))   # sin
s2 <- c(s2=function(x) esn(x, p=0.5))
s3 <- c(s3=function(x) esn(x, p=2.0))

sg <- function(x) psd(1/(1 + exp(-std(tcrossprod(x)))))


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

st <- function(x) dc(x, 'st', 0.05)
bn <- function(x) dc(x, 'bn')
ps <- function(x) dc(x, 'ps')
ex <- function(x) dc(x, 'ex')
ca <- function(x) dc(x, 'ca', 0.05)
ch <- function(x) dc(x, 'ch')
i1 <- function(x) x
i2 <- function(x) x^2 + x
i3 <- function(x) scale((x + 1)^3, scale=FALSE)
i2 <- function(x) scale((x + 1)^2, scale=FALSE)
i3 <- function(x) scale((x + 1)^3, scale=FALSE)
sn <- function(x) sin(2 * pi * x)


## genomic models
a1 <- ~ a
a2 <- ~ a + I(a^2)
aa <- ~ a + I(a^2) + a:a[]
ax <- ~ a:a[]

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
#' @param n the number of components
#' @param m the number of instances
#' @param s scale the components by this much.
get.vcs <- function(n=1, v=1)
{
    v <- rep(v, l=n)
    rchisq(n, v)
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
get.sim <- function(dat, frq=1, lnk=i1, eps=1, vcs=NULL, oks=p1, yks=oks, ...)
{
    dot <- list(...)
    mdl <- dot$mdl %||% a1
    svc <- dot$svc %||% 1
    nvc <- length(oks)
    vcs <- dot$vcs %||% get.vcs(nvc, svc)

    Q <- with(dat, length(unique(dvp[dvp > 0])))
    R <- with(dat, length(unique(evl[evl > 0])))
    one <- function(gmx)
    {
        P <- with(dat, NCOL(gmx))
        N <- with(dat, NROW(gmx))
        ## functional SNP mask
        fmk <- sample(c(rep(TRUE, P * frq), rep(FALSE, P - P * frq)))
        within(list(gmx=gmx),
        {
            fmx <- gsm(mdl, gmx[, fmk], rm.dup=FALSE)
            fnl <- krn(fmx, oks)
            eta <- drop(mvrnorm(1, rep(0, nrow(gmx)), cmb(fnl, vcs)[[1]]))
            rsp <- lnk(eta) + rnorm(nrow(gmx), 0, sqrt(eps))
            knl <- krn(gmx, yks)
        })
    }
    ## core effect
    dat <- within(dat, dvp <- one(gmx[dvp > 0, ]))
    dat <- within(dat, evl <- one(gmx[evl > 0, ]))
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
