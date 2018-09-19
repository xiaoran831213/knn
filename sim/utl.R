## kernels combinations
library(MASS)
source('R/kpl.R')
id <- c(id=function(x) diag(1, NROW(x)))
p1 <- c(o1=function(x) ply(x, degree=1))
p2 <- c(o1=function(x) ply(x, degree=1),
        o2=function(x) ply(x, degree=2))
p3 <- c(o1=function(x) ply(x, degree=1),
        o2=function(x) ply(x, degree=2),
        o3=function(x) ply(x, degree=3))
ga <- c(ga=function(x) gau(x))
lp <- c(lp=function(x) lap(x))
ib <- c(ib=function(x) ibs(x))

gk <- c(ga=function(x) gau(x),
        ki=function(x) kin(x))
kg <- c(ki=function(x) kin(x),
        ga=function(x) gau(x))

ki <- c(ki=function(x) kin(x))          # kinship
s1 <- c(s1=function(x) esn(x, p=1.0))   # sin
s2 <- c(s2=function(x) esn(x, p=0.5))
s3 <- c(s3=function(x) esn(x, p=2.0))

ks <- c(ki=function(x) kin(x),
        sn=function(x) esn(x, p=1.0))
sk <- c(sn=function(x) esn(x, p=1.0),
        ki=function(x) kin(x))

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
get.gmx <- function(gmx, N=NULL, P=NULL, Q=4, R=1)
{
    N <- N %||% as.integer(nrow(gmx) * .2)
    P <- min(P, ncol(gmx))

    dvp.msk <- c(rep(1L:Q, each=N), rep(0L, N * R))
    evl.msk <- c(rep(0L, N * Q), rep(1L:R, each=N))
    
    ## indices
    idx <- sample.int(nrow(gmx))[seq.int(N * Q + N * R)]
    jdx <- seq(sample.int(ncol(gmx) - P, 1L), l=P)
    gmx <- gmx[idx, jdx]

    ## remove degeneracy
    af <- colMeans(gmx) / 2.0
    gmx <- gmx[, pmin(af, 1 - af) >= 0.05]
    
    ## divide
    list(gmx=gmx, dvp=dvp.msk, evl=evl.msk)
}

#' get variance components randomly
#'
#' @param n the number of components
#' @param m drawn the components from rchisq or softmax
#' softmax ensure that the compoents adds up to one.
#' @param s scale the components by this much.
get.vcs <- function(n, s=rep(1, n), ...)
{
    s <- rep(s, l=n)
    e <- rchisq(n, s)
    e
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
#' corresponding, generted response.
get.sim <- function(dat, frq=1, lnk=i1, eps=1, vcs=NULL, oks=p1, yks=oks, ...)
{
    dot <- list(...)
    mdl <- dot$mdl %||% a1
    svc <- dot$svc %||% 1

    Q <- with(dat, length(unique(dvp[dvp > 0])))
    R <- with(dat, length(unique(evl[evl > 0])))
    P <- with(dat, NCOL(gmx))
    N <- with(dat, NROW(gmx))
    
    ## functional SNP mask
    fmk <- sample(c(rep(TRUE, P * frq), rep(FALSE, P - P * frq)))
    
    ## true variance components linking x to y, in log scale
    nvc <- length(oks)
    if(is.null(vcs))
        vcs <- get.vcs(nvc, 'r', svc)

    do.sim <- function(gmx)
    {
        within(list(gmx=gmx),
        {
            fmx <- gsm(mdl, gmx[, fmk], rm.dup=FALSE)
            fnl <- krn(fmx, oks)
            fcv <- cmb(fnl, vcs)[[1]]
            eta <- drop(mvrnorm(1, rep(0, nrow(gmx)), fcv, empirical=FALSE))
            rsp <- lnk(eta) + rnorm(nrow(gmx), 0, sqrt(eps))
            knl <- krn(gmx, yks)
        })
    }
    ## core effect
    dat <- within(dat, dvp <- do.sim(gmx[dvp > 0, ]))
    dat <- within(dat, evl <- do.sim(gmx[evl > 0, ]))
    dat
}

## randomly pick a RDS file from a directory
get.rds <- function(.)
{
    fs <- dir(., "[.]rds$", ful=TRUE)
    sample(fs, 1)
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

knl.lm <- function(knl, coef=NULL)
{
    N <- nrow(knl[[1]])
    knl <- c(list(diag(N)), knl)

    udx <- upper.tri(knl[[1]], TRUE)
    kdx <- seq(3, l=(length(knl) - 2))

    dat <- do.call(cbind, lapply(knl, `[`, udx))

    print(cor(dat))
    for(i in kdx)
    {
        m <- lm(dat[, i] ~ dat[, 1:(i-1)])
        dat[, i] <- m$residuals
    }

    print(cor(dat))
    for(i in kdx)
    {
        k <- knl[[i]]
        k[udx] <- dat[, i]
        k <- t(k)
        k[udx] <- dat[, i]
        knl[[i]] <- k
    }
    
    for(i in kdx)
    {
        knl[[i]] <- std(knl[[i]])
    }
    knl[-1]
}

knl.pc <- function(knl)
{
    udx <- upper.tri(knl[[1]], TRUE)
    dat <- do.call(cbind, lapply(knl, `[`, udx))
    pcs <- with(svd(dat), u %*% diag(d))

    for(i in seq_along(knl))
    {
        k <- knl[[i]]
        k[udx] <- pcs[, i]
        k <- t(k)
        k[udx] <- pcs[, i]
        knl[[i]] <- k
    }
    knl
}
upt <- function(x, ...) x[upper.tri(x), ...]
