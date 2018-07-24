## model aggregation

ccv <- function(rpt, ...)
{
    Q <- length(rpt)
    ds <- sprintf('d%02d', 1L:Q)
    ms <- sprintf('m%02d', 1L:Q)
    ret <- array(list(), c(Q, Q), list(ds, ms))
    err <- c()
    for(i in seq.int(Q))
    {
        ki <- rpt[[i]]$knl              # data[i]
        ri <- rpt[[i]]$rsp              # data[i]
        for(j in seq.int(Q))
        {
            pj <- rpt[[j]]$par          # model[j]
            ij <- vpd(ri, ki, pj, rt=0, ...)
            ret[i, j] <- list(ij)
        }

        ## collect the types of error/loss
        err <- union(err, names(ij))
    }

    ## organize by error type
    ret <- sapply(err, function(e)
    {
        apply(ret, c(1, 2), function(r) r[[1]][e])
    },
    simplify=FALSE)
    ret
}

cwt <- function(rpt, type='nlk', ...)
{
    ## relative generlization in training sets
    if(type == 'nlk')
    {
        tmp <- -sweep(rpt, 2L, diag(rpt))
        wmt <- exp(colSums(tmp))
        ## wmt <- exp(colSums(diag(rpt) - rpt))
    }
    if(type == 'cyh')
    {
        wmt <- exp(colSums(rpt - diag(rpt)))
    }
    if(type %in% c('mse', 'loo'))
    {
        tmp <- -sweep(rpt, 2L, diag(rpt))
        wmt <- exp(colSums(tmp))
    }
    if(type == 'ssz')
    {
        wmt <- exp(colSums(diag(rpt) - rpt))
    }
    wmt <- wmt / sum(wmt)
    wmt
}


#' aggregate cohorts
mat <- function(rpt, ...)
{
    ## report of leave one (cohort) out
    ## row: cohort, col: errors
    v <- ccv(rpt)
    w <- sapply(names(v), function(e)
    {
        cwt(v[[e]], e)
    })
    
    ## row: parameter, col: cohort
    p <- sapply(rpt, function(.) .$par)
    p %*% w
}
