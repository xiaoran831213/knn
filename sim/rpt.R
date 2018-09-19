library(dplyr)

sel <- dplyr::select
gba <-  function(x, ...)
{
    group_by_at(.tbl=x, .vars=vars(...))
}

## clean up parentheses
.cp <- function(x)
{
    x <- gsub("^c", "", x)
    x <- gsub("[ ()]", "", x)
    x <- gsub(",", "+", x)
    x
}

readRtm <- function(..., ref=TRUE)
{
    rpt <- list()
    dot <- c(...)
    for(f in dir(dot, '.rds', full=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next
        r <- subset(r, dat=='dvp' & key=='rtm', -dat)
        rpt <- c(rpt, list(r))
    }
    rpt <- do.call(rbind, rpt)
    
    if(!('lnk' %in% names(rpt)))
        rpt$lnk <- 'I'
    if(!('mdl' %in% names(rpt)))
        rpt$mdl <- 'a'
    if(!('oks' %in% names(rpt)))
        rpt$oks <- 'p'
    rpt <- within(rpt,
    {
        sim <- sprintf("%s(%s)~%s, x~%s", lnk, .cp(oks), .cp(yks), mdl)
    })

    ## the longest time span as measuring stick
    if(ref)
    {
        rpt <- as_tibble(rpt) %>% group_by(seed) %>% mutate(val=val/max(val)) %>% ungroup
    }
    rpt <- rpt %>% sel(-c(seed, oks, lnk, yks, mdl))
    as.data.frame(rpt)
}
    
readRpt <- function(..., ref=TRUE)
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
        if(ref)
        {
            ## MSE: null model as a measuring stick
            mse <- subset(evl, key == 'mse' & mtd != 'nul')
            nul <- subset(evl, key == 'mse' & mtd == 'nul', val)
            mse <- within(mse, val <- val / unlist(nul))

            ## NLK: null model as a measuring stick
            nlk <- subset(evl, key == 'nlk' & mtd != 'nul')
            nul <- subset(evl, key == 'nlk' & mtd == 'nul', val)
            nlk <- within(nlk, val <- exp(2 * (val - unlist(nul))))

            ## NLK: null model as a measuring stick
            loo <- subset(evl, key == 'loo' & mtd != 'nul')
            nul <- subset(evl, key == 'loo' & mtd == 'nul', val)
            loo <- within(loo, val <- val / unlist(nul))
            
            ## CYH: null model is meaning less
            cyh <- subset(evl, key == 'cyh' & mtd != 'nul')

            rpt <- c(rpt, list(rbind(mse, nlk, cyh, loo)))
        }
        else
        {
            rpt <- c(rpt, list(evl))
        }
    }
    rpt <- do.call(rbind, rpt)
    rpt <- subset(rpt, !is.infinite(val))
    
    ## corr(y, y_hat) for null model is meaning less
    rpt <- subset(rpt, !(key == 'cyh' & mtd == 'nul'))
    
    if(!('lnk' %in% names(rpt)))
        rpt$lnk <- 'I'
    if(!('mdl' %in% names(rpt)))
        rpt$mdl <- 'a'
    if(!('oks' %in% names(rpt)))
        rpt$oks <- 'p'
    rpt <- within(rpt,
    {
        sim <- sprintf("%s(%s)~%s, x~%s", lnk, .cp(oks), .cp(yks), mdl)
    })
    rpt <- subset(rpt, se=-c(seed, oks, lnk, yks, mdl))
    rpt
}

readBia <- function(..., ref=TRUE)
{
    rpt <- list()
    dot <- c(...)
    for(f in dir(dot, '.rds$', full=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next

        dvp <- subset(r, dat == 'dvp', -dat)
        rpt <- c(rpt, list(dvp))
    }
    rpt <- do.call(rbind, rpt)
    rpt <- subset(rpt, !is.infinite(val))
    
    if(!('lnk' %in% names(rpt)))
        rpt$lnk <- 'I'
    if(!('mdl' %in% names(rpt)))
        rpt$mdl <- 'a'
    if(!('oks' %in% names(rpt)))
        rpt$oks <- 'p'
    rpt <- within(rpt,
    {
        sim <- sprintf("%s(%s)~%s, x~%s", lnk, .cp(oks), .cp(yks), mdl)
    })
    rpt <- subset(rpt, se=-c(seed, oks, lnk, yks, mdl))
    rpt
}

get.unique <- function(d)
{
    i <- sapply(d, function(x) length(unique(x)) == 1)
    u <- as.list(d[1, i, drop=FALSE])
    t <- paste0(names(u), '=', u, collapse=', ')
    d <- d[, -which(i)]
    list(unq=u, ttl=t, dat=d)
}

## infer title from repeatative columns.
ttl <- function(d)
{
    i <- sapply(d, function(x) length(unique(x)) == 1)
    u <- as.list(d[1, i, drop=FALSE])
    paste0(names(u), '=', u, collapse=', ')
}

bxp <- function(dat, axi=val~mtd, ...)
{
    ..s <- list(...)
    for(. in names(..s)) assign(., ..s[[.]])

    ## labels
    vas <- all.vars(axi)
    ylab.bxp <- get0('ylab.bxp', inherits=FALSE, ifnotfound=vas[1])
    xlab.bxp <- get0('xlab.bxp', inherits=FALSE, ifnotfound=vas[2])

    ## limites
    ylim.bxp <- get0('ylim.bxp', inherits=FALSE, ifnotfound=c(0, 1))

    ## title
    title.bxp <- get0('title.bxp', inherits=FALSE, ifnotfound=get.unique(dat)$ttl)

    ## box plot
    boxplot(axi, dat, xlab=xlab.bxp, ylab=ylab.bxp, ylim=ylim.bxp, lex.order=FALSE, ...)
    abline(1, 0, col='red')
    title(title.bxp)
}


plt <- function(dat, oxy, ipt=NULL, out=NULL, ...)
{
    ..s <- list(...)

    ## main title
    .un <- get.unique(dat)
    ttl <- .un$ttl
    dat <- .un$dat

    ## inner plot
    if(is.null(ipt))
        ipt <- function(...) plot(NULL, NULL, 'n')
    if(is.function(ipt))
        ipt <- list(ipt)
    
    ## labels
    vas <- all.vars(oxy)
    olx <- get0('olx', inherits=FALSE, ifnotfound=vas[2])
    oly <- get0('oly', inherits=FALSE, ifnotfound=vas[1])

    ## groups
    xk <- if(olx == '.') NULL else as.factor(dat[, olx]) # x key
    yk <- if(oly == '.') NULL else as.factor(dat[, oly]) # y key
    xl <- if(is.null(xk)) NULL else unique(xk)
    yl <- if(is.null(yk)) NULL else unique(yk)
    nc <- if(is.null(xl)) 1 else length(xl)
    nr <- if(is.null(yl)) 1 else length(yl)

    ## prepair functions
    ipt <- rep(ipt, l=nc)

    gf <- list(x=xk, y=yk)
    gf <- gf[!sapply(gf, is.null)]
    gp <- split(dat, gf)
    
    ## plot
    ut <- 9 / max(nc, nr)
    if(!is.null(out))
        png(out, width=nc * ut, height=nr * ut, units='in', res=300)

    ## 
    par(mfrow=c(nr, nc), mar=c(2.25, 2.35, 1.2, 0), mgp=c(1.3, .5, 0))
    par(oma=c(0, 0, 2, 0) +.1)
    gc <- 1
    for(i in seq.int(l=nr))
    {
        for(j in seq.int(l=nc))
        {
            ipt[[j]](gp[[gc]], ...)
            gc <- gc + 1
        }
    }
    title(ttl, outer=TRUE, adj=0)

    if(!is.null(out))
        dev.off()
}
