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

        evl <- subset(r, dat == 'evl')
        dvp <- subset(r, dat == 'dvp')
        bia <- subset(r, dat == 'bia')
        
        ## MSE
        mse <- subset(evl, key == 'mse')
        nul <- subset(mse, mtd == 'NUL')$val
        mse <- within(mse, val <- val / nul)

        ## likelihood
        nlk <- subset(evl, key == 'nlk')
        nul <- subset(nlk, mtd == 'NUL')$val
        ldt <- within(subset(evl, key == 'ldt'), val <- exp(val - nul))
        yay <- within(subset(evl, key == 'yay'), val <- exp(val - nul))
        nlk <- within(nlk, val <- exp(val - nul))

        ## R^2
        rsq <- subset(evl, key=='rsq')

        ## running time
        rtm <- subset(dvp, key == 'rtm')
        rtm <- within(rtm, val <- val / median(val))
        rpt <- c(rpt, list(rbind(bia, rtm, mse, yay, ldt, nlk, rsq)))
    }
    rpt <- do.call(rbind, rpt)
    
    ## clean up parentheses
    .cp <- function(x) gsub("[~ ()]", "", x)
    rpt <- within(rpt,
    {
        oks <- .cp(oks)
        if(exists('lnk', inherits=FALSE))
        {
            oks <- sprintf("%s(%s)", lnk, sub('^~', '', oks))
            rm(lnk)
        }
        if(exists('yks', inherits=FALSE))
            yks <- .cp(yks)
        else
            yks <- "*"
        sim <- sprintf("%s~%s", oks, yks)
        
        if(exists('pss', inherits=FALSE) && length(unique(pss)) > 1)
            sim <- sprintf("%s, ps=%s", sim, pss)
        if(exists('bpt', inherits=FALSE) && length(unique(bpt)) > 1)
            sim <- sprintf("%s, bp=%s", sim, bpt)
        if(exists('wep', inherits=FALSE) && length(unique(wep)) > 1)
            sim <- sprintf("%s, ep=%s", sim, wep)
        if(exists('bmq', inherits=FALSE) && length(unique(bmq)) > 1)
            sim <- sprintf("%s, bm=%s", sim, bmq)
        if(length(unique(Q)) > 1)
            sim <- sprintf("%s, Q=%s", sim, Q)
    })
    ## rpt <- within(rpt, {H <- N * R; N <- N * Q})
    rpt <- within(rpt, rm(R, Q, oks, seed))
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
    d <- filter(d, !mtd %in% c('NUL'))
    r <- list()
    
    ## report: errors
    if(errs)
    {
        rpt <- filter(d, key %in% c('rsq', 'mse', 'nlk', 'ldt', 'yay'))
        rpt <- mutate(rpt, val=pmin(val, 1.99), val=pmax(val, 0.01))
        g <- ggplot(rpt, aes(x=mtd, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(sim~key)
        g <- g + ylim(c(0, 2))
        g <- g + ggtitle(ttl(rpt))
        f <- paste0(file.path('~/img', name), '.rpt.png'); print(f)
        ggsave(f, g, width=16, height=10)
        r <- c(r, rpt=g)
    }

    if(bias)
    {
        ## report: bias
        bia <- filter(d, dat=='bia')
        bia <- mutate(bia, key=sub('^.*[.]', '', key))
        bia <- mutate(bia, val=pmin(val, 1.99), val=pmax(val, -1.99))
        g <- ggplot(bia, aes(x=key, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(sim~mtd)
        g <- g + ylim(c(-2, 2))
        g <- g + ggtitle(ttl(bia))
        f <- paste0(file.path('~/img', name), '.bia.png'); print(f)
        ggsave(f, g, width=16, height=10)
        r <- c(r, bia=g)
    }

    ## return
    invisible(r)
}

rpt1 <- function(cache=FALSE)
{
    ## prpt('b00', cache, bias=TRUE)
    ## prpt('b01', cache, bias=TRUE)
    ## prpt('b02', cache, bias=TRUE)
    ## prpt('b03', cache, bias=TRUE)
    ## prpt('b10', cache, bias=TRUE)
    ## prpt('b11', cache, bias=TRUE)
    ## prpt('b15', cache, bias=TRUE)
    ## prpt('b16', cache, bias=TRUE)
    ## prpt('b20', cache, bias=TRUE)
    ## prpt('b21', cache, bias=TRUE)
    ## prpt('b25', cache, bias=TRUE)
    ## prpt('b26', cache, bias=TRUE)
}

rpt2 <- function(cache=FALSE, bias=TRUE)
{
    prpt('a00', cache, bias=bias)
    prpt('a01', cache, bias=bias)
}

rpt3 <- function(cache=FALSE, bias=TRUE)
{
    prpt('d00', cache, bias=bias)
    prpt('d01', cache, bias=bias)
}
