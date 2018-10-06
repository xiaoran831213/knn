library(ggplot2)
library(dplyr)

plt.cor <- function(dat)
{
    d1 <- list()
    nm <- names(dat)
    for(i in seq_along(dat))
    {
        for(j in seq_along(dat))
        {
            p1 <- data.frame(dx=dat[, i], dy=dat[, j], vx=nm[i], vy=nm[j])
            d1 <- c(d1, list(p1))
        }
    }
    d1 <- do.call(rbind, d1)

    gc <- ggplot(d1, aes(x=dx, y=dy))
    gc <- gc + geom_point(size=.5)
    gc <- gc + facet_grid(vx~vy)

    th <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    gc <- gc + th
    gc
}

S01 <- function(r, k, ...) within(subset(r, key==k), {val <- val - min(val); val <- val / max(val)})
S02 <- function(r, k, ...)
{
    r <- subset(r, key == k)
    within(r, val <- val / subset(r, mtd=='nul')$val)
}
S03 <- function(r, k, ...) within(subset(r, key == k), val <- val / val[1])

sel <- dplyr::select

readRPT <- function(..., ref=1)
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

        mse <- S02(evl, 'mse')
        nlk <- S02(evl, 'nlk')
        cyh <- within(subset(evl, key=='cyh'), val <- val^2)
        rtm <- S03(dvp, 'rtm')
        rpt <- c(rpt, list(rbind(bia, rtm, mse, nlk, cyh)))
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
        dat <- readRPT(fdr, ref=2)
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
        rpt <- mutate(rpt, val=pmin(val, 1.99))
        g <- ggplot(rpt, aes(x=mtd, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(sim~key)
        g <- g + ylim(c(0, 2))
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

rpt1 <- function(cache=FALSE)
{
    prpt('d00', cache, bias=FALSE)
    prpt('d01', cache, bias=FALSE)

    prpt('h00', cache, bias=FALSE)
    prpt('h01', cache, bias=FALSE)
    prpt('h02', cache, bias=FALSE)
    prpt('h03', cache, bias=FALSE)

    prpt('h10', cache, bias=FALSE)
    prpt('h11', cache, bias=FALSE)
    prpt('h12', cache, bias=FALSE)
    prpt('h13', cache, bias=FALSE)
}

rpt2 <- function(cache=FALSE)
{
    prpt('a00', cache, bias=FALSE)
    prpt('a01', cache, bias=FALSE)
    prpt('a02', cache, bias=FALSE)
    prpt('a03', cache, bias=FALSE)

    prpt('a10', cache, bias=FALSE)
    prpt('a11', cache, bias=FALSE)
    prpt('a12', cache, bias=FALSE)
    prpt('a13', cache, bias=FALSE)
}

rpt3 <- function(cache=FALSE)
{
    prpt('p00', cache, bias=FALSE)
    prpt('p01', cache, bias=FALSE)
    prpt('p02', cache, bias=FALSE)
    prpt('p03', cache, bias=FALSE)

    prpt('p10', cache, bias=FALSE)
    prpt('p11', cache, bias=FALSE)
    prpt('p12', cache, bias=FALSE)
    prpt('p13', cache, bias=FALSE)
}

