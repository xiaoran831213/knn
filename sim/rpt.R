library(ggplot2)
library(dplyr)

S01 <- function(r, k, ...)
{
    ## r <- subset(r, key==k & mtd!='nul' & mtd!='fun')
    r <- subset(r, key==k)
    r <- within(r,
    {
        val <- val - min(val)
        val <- val / max(val)
    })
    r
}

sel <- dplyr::select

readRPT <- function(..., ref=TRUE)
{
    rpt <- list()
    dot <- c(...)
    for(f in dir(dot, '.rds$', full=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next

        evl <- subset(r, dat == 'evl', -dat)
        dvp <- subset(r, dat == 'dvp', -dat)
        bia <- subset(dvp, grepl('^bia[.]', key))
        if(ref)
        {
            mse <- S01(evl, 'mse')
            nlk <- S01(evl, 'nlk')
            ## cyh <- S01(evl, 'cyh')
            cyh <- subset(evl, key=='cyh')
            cyh <- within(cyh, val <- val^2)
            rtm <- S01(dvp, 'rtm')
            rpt <- c(rpt, list(rbind(bia, rtm, mse, nlk, cyh)))
        }
        else
        {
            rtm <- subset(dvp, key=='rtm')
            rpt <- c(rpt, list(rbind(bia, rtm, evl)))
        }
    }
    rpt <- do.call(rbind, rpt)
    rpt <- subset(rpt, !is.infinite(val))

    ## clean up parentheses
    .cp <- function(x) gsub("[~ ()]", "", x)
    rpt <- within(rpt,
    {
        oks <- .cp(oks)
        if(exists('lnk', inherits=FALSE))
            oks <- sprintf("%s(%s)", lnk, sub('^~', '', oks))
        sim <- sprintf("%s~%s", oks, .cp(yks))
    })
    rpt <- within(rpt, {H <- N * R; N <- N * Q})
    rpt <- within(rpt, rm(R, Q, oks, yks, seed))
    rpt
}

## infer title from repeatative columns.
ttl <- function(d)
{
    i <- sapply(d, function(x) length(unique(x)) == 1)
    u <- as.list(d[1, i, drop=FALSE])
    paste0(names(u), '=', u, collapse=', ')
}

## fetch reports
d0 <- function(name, use.cache=TRUE)
{
    fdr <- file.path('sim/run', name)
    rds <- paste0(fdr, '.rds')
    if(file.exists(rds) && use.cache)
        dat <- readRDS(rds)
    else
    {
        dat <- readRPT(fdr, ref=TRUE)
        saveRDS(dat, rds)
    }
    print(rds)
    dat
}

## plot reports
prpt <- function(name, cache=TRUE, errs=TRUE, bias=TRUE)
{
    d <- as_tibble(d0(name, cache))
    d <- filter(d, !mtd %in% c('fun', 'nul'))
    r <- list()

    ## report: errors
    if(errs)
    {
        rpt <- filter(d, key %in% c('cyh', 'rtm', 'mse', 'nlk'))
        g <- ggplot(rpt, aes(x=mtd, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(sim~key)
        g <- g + ggtitle(ttl(rpt))
        f <- paste0(file.path('~/img', name), '.rpt.png'); print(f)
        ggsave(f, g, width=10, height=11)
        r <- c(r, rpt=g)
    }

    if(bias)
    {
        ## report: bias
        bia <- filter(d, grepl('^bia[.]', key))
        bia <- mutate(bia, key=sub('^.*[.]', '', key))
        g <- ggplot(bia, aes(x=key, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(sim~mtd)
        g <- g + ggtitle(ttl(bia))
        f <- paste0(file.path('~/img', name), '.bia.png'); print(f)
        ggsave(f, g, width=10, height=11)
        r <- c(r, bia=g)
    }

    ## return
    invisible(r)
}

rpt1 <- function()
{
    prpt('d00', FALSE, bias=FALSE)
    prpt('d01', FALSE, bias=FALSE)
    prpt('d02', FALSE, bias=FALSE)
    prpt('d03', FALSE, bias=FALSE)
}
