#' get a consecutive section from the genomic matrix
#' for training and testing dataset.
#' 
#' @param the basic gnomic matrix
#' @param N the size of development data
#' @param H the size of evaluation data
#' @param P the number of variants
get.gmx <- function(gmx, N, H=N, P=ncol(gmx) * .5)
{
    H <- min(N, nrow(gmx) - N)
    P <- min(P, ncol(gmx))
    
    ## individual indices
    idx <- sample.int(nrow(gmx), N + H)
    ## variant indices
    jdx <- seq(sample.int(ncol(gmx) - P, 1), l=P)

    ## developing data
    gmx.dvp <- as.matrix(gmx[idx[+(1:N)], jdx])

    ## evaluating data
    gmx.evl <- as.matrix(gmx[idx[-(1:N)], jdx])
    
    list(dvp=gmx.dvp, evl=gmx.evl)
}

sample.gmx <- function(gmx, N=NULL, P=NULL, Q=3)
{
    N <- N %||% as.integer(nrow(gmx) * .2)
    P <- P %||% as.integer(ncol(gmx) * .5)

    ## indices
    idx <- sample.int(nrow(gmx))[seq.int(N * Q)]
    jdx <- seq(sample.int(ncol(gmx) - P, 1), l=P)
    
    ret <- lapply(1:Q, function(i)
    {
        as.matrix(gmx[seq.int(N * i, l=N), jdx])
    })
    ret
}

sample.vcs <- function(e, k, rep=3)
{
    replicate(rep, c(e, rchisq(k - 1, 1)), FALSE)
}

get.sim <- function(gmx, vcs, frq=1, lnk=I, oks=c(id, p1), ejt=.1)
{
    P <- NCOL(gmx)
    N <- NROW(gmx)
    k <- length(oks)

    ## functional SNP mask
    fmk <- sample(c(rep(1, P * frq), rep(0, P - P * frq)))

    ## true variance components linking x to y, in log scale
    vcs <- exp(log(vcs) + rnorm(k, 0, ejt))
    eps <- vcs[1]

    ## generate
    fmx <- sweep(gmx, 2L, fmk, `*`)
    fnl <- krn(fmx, oks)
    knl <- krn(gmx, oks)
    fcv <- cmb(fnl, vcs)[[1]]
    kcv <- cmb(knl, vcs)[[1]]
    rsp <- mvn(1, fcv) %>% drop %>% lnk

    ## oracle fit and null fit:
    rpt <- list()
    rpt <- cl(rpt, DF(mtd='fmx', lmm(rsp, fcv, eps)))
    rpt <- cl(rpt, DF(mtd='gmx', lmm(rsp, kcv, eps)))
    rpt <- cl(rpt, DF(mtd='nul', nul(rsp)))
    rpt <- do.call(rbind, rpt)

    ## return
    list(gmx=gmx, rpt=rpt, rsp=rsp, fcv=fcv, kcv=kcv, vcs=vcs, eps=eps)
}
