## simulation utilities
source('R/utl.R')

pars <- function(dvp, ref=NULL)
{
    par <- dvp %$% 'par'
    key <- unlist(lapply(par, names), use.names=FALSE)
    if(!is.null(ref))
        key <- c(names(ref), key)
    key <- unique(key)
    mtd <- names(dvp)

    par <- sapply(par, `[`, key)
    rownames(par) <- key
    par[is.na(par)] <- 0
    par <- data.frame(t(par))
    
    if(!is.null(ref))
    {
        ref <- ref[key]
        names(ref) <- key
        ref[is.na(ref)] <- 0
        par <- rbind(par, REF=ref)
    }
    par
}

bias <- function(dvp, eps, vcs)
{
    ## bias assesment
    par <- dvp %$% 'par'                # estimates
    ref <- c(EPS=eps, vcs)              # reference

    ## parameter names, and methods
    key <- unique(unlist(c(names(ref), sapply(par, names)), use.names=FALSE))
    mtd <- names(dvp)
    
    ## alignment:
    ## 1) the estimates
    par <- sapply(par, `[`, key)
    rownames(par) <- key
    par[is.na(par)] <- 0

    ## 2) reference
    ref <- ref[key]
    names(ref) <- key
    ref[is.na(ref)] <- 0

    ## the table of bias
    bia <- data.frame(t(par - ref))
    bia <- reshape(bia, key, 'val', timevar='key', idvar='mtd', ids=mtd, times=key, direction='l')
    DF(dat='bia', bia)
}
