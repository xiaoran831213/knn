#' get a consecutive section from the genomic matrix
#' for training and testing dataset.
#' 
#' @param gmx the gnomic matrix
#' @param N the size of development data
#' @param P the number of variants
get.gmx <- function(gmx, N=NULL, P=NULL, Q=4, R=1)
{
    N <- N %||% as.integer(nrow(gmx) * .2)
    P <- min(P, ncol(gmx))
    
    ## indices
    idx <- sample.int(nrow(gmx))[seq.int(N * Q * 2)]
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

get.vcs <- function(n, mtd=c('softmax', 'rchisq'), sc=1)
{
    mtd <- match.arg(mtd, c('softmax', 'rchisq'))
    if(mtd == 'softmax')
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

get.sim <- function(G, frq=1, lnk=I, eps=1, V=p1, het=.1, vc=NULL)
{
    if(!is.list(G))
        G <- list(G)

    P <- NCOL(G[[1]])
    N <- NROW(G[[1]])
    k <- length(V)
    
    ## functional SNP mask
    fmk <- sample(c(rep(1, P * frq), rep(0, P - P * frq)))
    
    ## true variance components linking x to y, in log scale
    nvc <- length(V)
    if(is.null(vc))
        vcs <- c(eps=eps, vc=get.vcs(nvc, 'r', 2))
    else
        vcs <- c(eps=eps, vc=vc)
    cvs <- c(id, V)

    ## jittering
    het <- rep(het, l=length(G))
    
    ## effects
    ## eft <- rnorm(P)

    ret <- mapply(function(gmx, a)
    {
        .vc <- vcs * (1 - a) + c(eps, get.vcs(nvc, 'r', 2)) * a

        ## generate
        ## fmx <- sweep(gmx, 2L, eft, `*`)
        fmx <- sweep(gmx, 2L, fmk, `*`)
        ## fmx <- gmx[, as.logical(fmk)]
        fnl <- krn(fmx, cvs)
        knl <- krn(gmx, cvs)
        fcv <- cmb(fnl, .vc)[[1]]
        kcv <- cmb(knl, .vc)[[1]]
        rsp <- mvn(1, fcv) %>% drop %>% lnk

        ## oracle fit and null fit:
        rpt <- list()
        rpt <- CL(rpt, DF(mtd='fmx', lmm(rsp, fcv, eps)))
        rpt <- CL(rpt, DF(mtd='gmx', lmm(rsp, kcv, eps)))
        rpt <- CL(rpt, DF(mtd='nul', nul(rsp)))
        rpt <- do.call(rbind, rpt)

        ## return
        list(gmx=gmx, rpt=rpt, rsp=rsp, vcs=.vc)
    },
    G, het, SIMPLIFY=FALSE)
    ret
}

get2 <- function(G, frq=1, lnk=I, eps=1, V=p1, het=.1, vc1=NULL)
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
        vc1 <- get.vcs(nvc, 'r', 2)
    vcs <- c(eps=eps, vc=vc1)
    cvs <- c(id, V)

    ## jittering
    het <- rep(het, l=length(G))
    
    jit <- lapply(G, function(gmx)
    {
        .vc <- c(eps, get.vcs(nvc, 'r', 2))

        ## generate
        eft <- rnorm(P)

        fmx <- sweep(gmx, 2L, eft, `*`)
        ## fmx <- sweep(gmx, 2L, fmk, `*`)
        fmx <- gmx[, as.logical(fmk)]
        fnl <- krn(fmx, cvs)
        fcv <- cmb(fnl, .vc)[[1]]

        mvn(1, fcv) %>% drop %>% lnk
    })
    whl <- with(list(),
    {
        gmx <- do.call(rbind, G)
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
        rsp <- whl[[i]] * (1-het[i]) + jit[[i]] * het[i]
        mix[[i]] <- list(rsp=rsp, gmx=G[[i]])
    }

    mix
}

