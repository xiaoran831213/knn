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

## link functions
#' Distribution Convertor
#'
#' covert an assumed empirical nromal distribution to another distribution.
#'
#' x: the data that descript the empirical normal
#' q: the quantile function of the target distribution
#' .: additional paramters required by the target quantile (e.g., degree
#' of freedom for t and chisq, min and min for unif, etc.)
dc <- function(x, d=NULL, ...)
{
    v <- var(x)                         # variance
    s <- sd(x)                          # sd
    m <- mean(x)                        # mean
    p <- pnorm(x, m, s)                 # percentile

    y <- switch(
        d,
        bn={2 * v},
        ps={v},
        ex={1 / s},
        ch={v},
        0)

    z <- switch(
        d,
        st={     qt(p, (2 * v) * 2/ max(v - 1, 0))},
        bn={ qbinom(p, ceiling(4 * v), .5)},
        ps={  qpois(p, v)},
        ex={   qexp(p, 1 / s)},
        ca={qcauchy(p) -> .; v * . * 8/ sd(.)},
        ch={ qchisq(p, v / 2)})

    cdc <- as.logical(get0('cdc', ifnotfound=FALSE))
    print(list(cdc=cdc))
    if(cdc)
        z <- z - y
    z
}

st <- function(x) dc(x, 'st')
bn <- function(x) dc(x, 'bn')
ps <- function(x) dc(x, 'ps')
ex <- function(x) dc(x, 'ex')
ca <- function(x) dc(x, 'ca')
ch <- function(x) dc(x, 'ch')
i1 <- function(x) x
i2 <- function(x) scale((x + 1)^2, scale=FALSE)
i3 <- function(x) scale((x + 1)^3, scale=FALSE)
i2 <- function(x) scale((x + 1)^2, scale=FALSE)
i3 <- function(x) scale((x + 1)^3, scale=FALSE)
sn <- function(x) sin(2 * pi * x)
sg <- function(x) 1/(1 + exp(-x))

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
    
    ## indices
    idx <- sample.int(nrow(gmx))[seq.int(N * Q + N * R)]
    jdx <- seq(sample.int(ncol(gmx) - P, 1), l=P)
    gmx <- gmx[idx, jdx]

    ## remove degeneracy
    af <- colMeans(gmx) / 2.0
    gmx <- gmx[, pmin(af, 1 - af) >= 0.05]
    
    ## divide
    dvp <- lapply(seq(0, l=Q), function(i)
    {
        as.matrix(gmx[seq.int(1 + N * i, l=N), ])
    })
    names(dvp) <- sprintf('d%02d', seq_along(dvp))

    evl <- lapply(seq(Q, l=R), function(i)
    {
        as.matrix(gmx[seq.int(1 + N * i, l=N), ])
    })
    names(evl) <- sprintf('e%02d', seq_along(evl))
    
    list(dvp=dvp, evl=evl)
}

#' get variance components randomly
#'
#' @param n the number of components
#' @param m drawn the components from rchisq or softmax
#' softmax ensure that the compoents adds up to one.
#' @param s scale the components by this much.
get.vcs <- function(n, m=c('rchisq', 'softmax'), svc=2, ...)
{
    m <- match.arg(m, c('softmax', 'rchisq'))
    if(m == 'softmax')
    {
        e <- rnorm(n)
        e <- exp(e)
        e <- e / sum(e) * svc
    }
    else
    {
        e <- rchisq(n, svc)
    }
    e
}

#' get simulated response
#'
#' @param G the a list of genotype data for different cohorts
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
get.sim <- function(G, frq=1, lnk=I, eps=1, V=p1, het=.1, svc=1, ...)
{
    dot <- list(...)
    mdl <- dot$mdl
    vc1 <- dot$vc1
    ## 1) het population
    if(!is.list(G))
        G <- list(G)

    Q <- length(G)
    P <- NCOL(G[[1]])
    N <- NROW(G[[1]])
    k <- length(V)
    
    ## functional SNP mask
    fmk <- sample(c(rep(1, P * frq), rep(0, P - P * frq)))
    
    ## true variance components linking x to y, in log scale
    nvc <- length(V)
    if(is.null(vc1))
        vc1 <- get.vcs(nvc, 'r', svc)
    vcs <- c(eps=0, vc=vc1)
    cvs <- c(id, V)

    ## jittering
    jit <- NULL
    if(any(het > 0))
    {
        het <- rep(het, l=length(G))
        jit <- lapply(G, function(gmx)
        {
            .vc <- get.vcs(nvc + 1, 'r', svc)

            ## generate
            eft <- rnorm(P)
            fmx <- sweep(gmx, 2L, eft, `*`)
            ## fmx <- sweep(gmx, 2L, fmk, `*`)
            fmx <- gmx[, as.logical(fmk)]
            fnl <- krn(fmx, cvs)
            fcv <- cmb(fnl, .vc)[[1]]

            rnorm(1) + mvn(1, fcv) %>% drop %>% lnk
        })
    }

    whl <- with(list(),
    {
        gmx <- do.call(rbind, G)
        ## fmx <- sweep(gmx, 2L, fmk, `*`)
        fmx <- gmx[, as.logical(fmk)]
        if('mdl' %in% names(dot))
            fmx <- gsm(G=fmx, ...)
        fnl <- krn(fmx, cvs)
        fcv <- cmb(fnl, vcs)[[1]]
        ## rsp <- mvn(1, fcv) %>% drop %>% lnk
        rsp <- mvrnorm(1, mu=rep(0, NROW(fcv)), Sigma=fcv) %>% drop %>% lnk
        rsp <- rsp + rnorm(NROW(rsp), 0, sqrt(eps))

        s1 <- cumsum(sapply(G, NROW))
        s0 <- c(0, s1[-Q]) + 1
        mapply(function(a, b)
        {
            rsp[a:b]
        },
        s0, s1, SIMPLIFY=FALSE)
    })
    
    mix <- list()
    for(i in seq.int(1, l=Q))
    {
        rsp <- whl[[i]]
        if(!is.null(jit))
            rsp <- rsp * (1 - het[i]) + jit[[i]] * het[i]
        mix[[i]] <- list(rsp=rsp, gmx=G[[i]])
    }
    mix
}

## randomly pick a RDS file from a directory
get.rds <- function(.)
{
    fs <- dir(., "[.]rds$", ful=TRUE)
    sample(fs, 1)
}
