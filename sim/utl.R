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

sample.gmx <- function(gmx, N=NULL, P=NULL, Q=4, H=N)
{
    N <- N %||% as.integer(nrow(gmx) * .2)
    P <- min(P, ncol(gmx))
    
    ## indices
    idx <- sample.int(nrow(gmx))[seq.int(N * Q + H)]
    jdx <- seq(sample.int(ncol(gmx) - P, 1), l=P)
    gmx <- gmx[idx, jdx]
    
    ## divide
    dvp <- lapply(seq(0, l=Q), function(i)
    {
        as.matrix(gmx[seq.int(1 + N * i, l=N), ])
    })
    names(dvp) <- sprintf('d%02d', seq_along(dvp))

    evl <- as.matrix(gmx[seq.int(N * Q + 1, l=H), ])
    
    c(dvp, list(evl=evl))
}

sample.vcs <- function(e, k, rep=3)
{
    replicate(rep, c(e, rchisq(k - 1, 1)), FALSE)
}

get.sim <- function(gms, vcs, frq=1, lnk=I, oks=c(p1), ejt=.1)
{
    if(!is.list(gms))
    {
        gms <- list(gms)
    }

    P <- NCOL(gms[[1]])
    N <- NROW(gms[[1]])
    k <- length(oks)

    ## functional SNP mask
    fmk <- sample(c(rep(1, P * frq), rep(0, P - P * frq)))
    
    ## true variance components linking x to y, in log scale
    ## print(vcs)
    nvc <- length(vcs)

    ## noise
    eps <- vcs[1]

    ## jittering
    ejt <- rep(ejt, l=length(gms))

    ret <- mapply(function(gmx, a)
    {
        ## jit <- rchisq(nvc - 1, 1)
        ## vcs[-1] <- vcs[-1] * (1 - a) + jit * a
        jit <- rchisq(nvc, 1)
        vcs <- vcs * (1 - a) + jit * a
    
        ## generate
        ## fmx <- sweep(gmx, 2L, fmk, `*`)
        fmx <- gmx[, as.logical(fmk)]
        fmx <- sweep(fmx, 2L, rnorm(sum(fmk)), `*`)
        fnl <- krn(fmx, oks)
        knl <- krn(gmx, oks)
        fcv <- cmb(fnl, vcs)[[1]]
        kcv <- cmb(knl, vcs)[[1]]
        rsp <- mvn(1, fcv) %>% drop %>% lnk
        rsp <- rsp - mean(rsp)

        ## oracle fit and null fit:
        rpt <- list()
        rpt <- cl(rpt, DF(mtd='fmx', lmm(rsp, fcv, eps)))
        rpt <- cl(rpt, DF(mtd='gmx', lmm(rsp, kcv, eps)))
        rpt <- cl(rpt, DF(mtd='nul', nul(rsp)))
        rpt <- do.call(rbind, rpt)

        ## return
        list(gmx=gmx, rpt=rpt, rsp=rsp, fcv=fcv, kcv=kcv, vcs=vcs, eps=eps)
    }, gms, ejt, SIMPLIFY=FALSE)
    ret
}
