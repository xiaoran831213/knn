## link functions

## genetic effect model
gsm <- function(mdl=~ a + d + r, G,   # model, and raw data.
                max.gvr=Inf,          # maximum number of variants
                max.tms=Inf,          # naximum number of terms
                rm.nic=TRUE,          # remove non-informative columns
                rm.dup=TRUE,          # remove duplicated columns
                ...)
{
    ## pick genomic variables
    if(ncol(G) > max.gvr)
        G <- G[, sort(sample.int(ncol(G), max.gvr))]

    ## pick out informative variants
    . <- apply(G, 2L, sd) > 0
    if(ncol(G) > sum(.))
    {
        cat('drop', ncol(G) - sum(.),'non-informative variants.\n')
        G <- G[., , drop=FALSE]
    }

    ## name the variants
    if(is.null(names(G)))
        colnames(G) <- sprintf('%02d', 1L:ncol(G))

    ## terms and varialbes
    tm <- terms(mdl)
    attr(tm, 'intercept') <- 0
    attr(tm, 'response') <- 0
    e <- all.vars(tm)
    if('g' %in% e)
        e <- c('a', 'd', 'r')
    
    ## genetic bases, additive, dominent, resessive
    b <- lapply(e, `==`, sample(rep(e, l=ncol(G))))
    names(b) <- e
    b <- lapply(b, function(.) G[, .])
    b$d <- 2 * (b$d > 0)
    b$r <- 2 * (b$r > 1)
    b$g <- with(b, cbind(a, d, r)[, colnames(G)])
    
    ## prepare model matrix
    m <- model.matrix(tm, b)

    ## fix term names
    n <- gsub('[I()]|\\[|\\]', '', colnames(m))
    n <- sapply(strsplit(n, ':'), function(x) paste0(sort(x), collapse=':'))
    n <- sub('([a-z]+)(\\^[0-9])([0-9]+)', '\\1\\3\\2', n)
    
    ## remove duplicated terms
    a <- attr(m, 'assign')
    ## . <- !duplicated(n)
    ## if(ncol(m) > sum(.))
    ## {
    ##     cat('drop', ncol(m) - sum(.), 'duplicated terms.\n')
    ##     m <- m[, ., drop=FALSE]
    ##     a <- a[.]
    ##     n <- n[.]
    ## }

    ## remove artificially created polynomial terms
    . <- !grepl('(.*):\\1', n)
    if(ncol(m) > sum(.))
    {
        cat('drop', ncol(m) - sum(.), 'artificial polynomial.\n')
        m <- m[, ., drop=FALSE]
        a <- a[.]
        n <- n[.]
    }

    ## assign decided term names
    colnames(m) <- n

    ## keep informative columns
    . <- apply(m, 2L, sd) > 0
    if(ncol(m) > sum(.) && rm.nic)
    {
        cat('drop', ncol(m) - sum(.), 'non-informative terms.\n')
        m <- m[, ., drop=FALSE]
        a <- a[.]
    }

    ## keep unique columns
    ## u <- unique(m, MARGIN=2L)
    ## if(ncol(m) > ncol(u))
    ## {
    ##     cat('drop', ncol(m) - ncol(u), 'duplicated columns.\n')
    ##     m <- u
    ## }

    ## crop excessive columns
    if(!is.null(max.tms) && ncol(m) > max.tms)
    {
        cat('drop', ncol(m) - max.tms, 'excessive columns.\n')
        . <- sort(sample.int(ncol(m), max.tms))
        m <- m[, ., drop=FALSE]
        a <- a[.]
    }
    attr(m, 'assign') <- a
    m
}


mdl.str <- function(mdl)
{
    mdl <- tail(as.character(mdl), 1)
    mdl <- gsub('[I() ]|\\[|\\]', '', mdl)
    mdl <- gsub('\\+-', '-', mdl)
    mdl
}

test.gsm <- function(N=9L, P=6L, ...)
{
    ## genotype
    G <- matrix(sample.int(3L, P * N, TRUE) - 1L, N, P)
    colnames(G) <- sprintf('%02d', 1L:P)

    ## model
    M <- ~ a

    ## expanded input matrix
    X <- gsm(M, G=G, rm.nic=TRUE, ...)

    ## return
    list(G=G, X=X, M=mdl.str(M))
}

