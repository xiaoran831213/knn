## kernels combinations
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
source('R/utl.R')
st <- function(x) ndc(x, qt, df=10)
bn <- function(x) ndc(x, qbinom, size=1, prob=.5)
ps <- function(x) ndc(x, qpois, lambda=1)
ex <- function(x) ndc(x, qexp, rate=1)
ca <- function(x) ndc(x, qcauchy, scale=.1)
ch <- function(x) ndc(x, qchisq, df=1)
i1 <- function(x) x
i2 <- function(x) scale((x + 1)^2)
i3 <- function(x) scale((x + 1)^3)
sn <- function(x) sin(2 * pi * x)
sg <- function(x) 1/(1 + exp(-x))

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
get.vcs <- function(n, m=c('rchisq', 'softmax'), sc=2, ...)
{
    m <- match.arg(m, c('softmax', 'rchisq'))
    if(m == 'softmax')
    {
        e <- rnorm(n)
        e <- exp(e)
        e <- e / sum(e) * sc
    }
    else
    {
        e <- rchisq(n, sc)
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
#' @param vc1 instead of randomly draw, use these variance compoent for the shared effect,
#' the heterogeneity effect is still randomly generated.
#'
#' @return a list of lists, where the inner lists contains original genotypes and
#' corresponding, generted response.
get.sim <- function(G, frq=1, lnk=I, eps=1, V=p1, het=.1, vc1=NULL, sc=1)
{
    ## 1) het population
    if(!is.list(G))
        G <- list(G)

    Q <- length(G)
    P <- NCOL(G[[1]])
    k <- length(V)
    
    ## functional SNP mask
    fmk <- sample(c(rep(1, P * frq), rep(0, P - P * frq)))
    
    ## true variance components linking x to y, in log scale
    nvc <- length(V)
    if(is.null(vc1))
        vc1 <- get.vcs(nvc, 'r', sc)
    vcs <- c(eps=eps, vc=vc1)
    cvs <- c(id, V)

    ## jittering
    if(any(het > 0))
    {
        het <- rep(het, l=length(G))
        jit <- lapply(G, function(gmx)
        {
            ## .vc <- c(eps, get.vcs(nvc, 'r', 2))
            .vc <- get.vcs(nvc + 1, 'r', sc)

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
    else
    {
        jit <- NULL
    }

    whl <- with(list(),
    {
        gmx <- do.call(rbind, G)
        ## fmx <- sweep(gmx, 2L, fmk, `*`)
        fmx <- gmx[, as.logical(fmk)]
        fnl <- krn(fmx, cvs)
        fcv <- cmb(fnl, vcs)[[1]]
        rsp <- mvn(1, fcv) %>% drop %>% lnk

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
