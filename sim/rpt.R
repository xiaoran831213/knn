library(dplyr)

## clean up parentheses
.cp <- function(x)
{
    x <- gsub("^c", "", x)
    x <- gsub("[ ()]", "", x)
    x <- gsub(",", "+", x)
    x
}

readRtm <- function(...)
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

    rpt <- within(rpt,
    {
        sim <- sprintf("%s(%s)~%s", lnk, .cp(oks), .cp(yks))
    })

    ## the longest time span as measuring stick
    rpt <- as_tibble(rpt) %>% group_by(seed) %>% mutate(val=val/max(val)) %>% ungroup
    rpt <- rpt %>% select(-c(seed, oks, lnk, yks))
    as.data.frame(rpt)
}
    
readRpt <- function(...)
{
    rpt <- list()
    dot <- c(...)
    for(f in dir(dot, '.rds$', full=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next
        
        ## use null model as a measuring stick
        evl <- subset(r, dat == 'evl', -dat)
        mse <- subset(evl, key == 'mse' & mtd != 'nul')
        nul <- subset(evl, key == 'mse' & mtd == 'nul', val)
        mse <- within(mse, val <- val / unlist(nul))

        nlk <- subset(evl, key == 'nlk' & mtd != 'nul')
        nul <- subset(evl, key == 'nlk' & mtd == 'nul', val)
        nlk <- within(nlk, val <- val / unlist(nul))

        rpt <- c(rpt, list(rbind(mse, nlk)))
    }
    rpt <- do.call(rbind, rpt)
    
    rpt <- within(rpt,
    {
        sim <- sprintf("%s(%s)~%s", lnk, .cp(oks), .cp(yks))
    })
    rpt <- subset(rpt, se=-c(seed, oks, lnk, yks))
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

bxp <- function(dat, axi=val~mtd, ...)
{
    ..s <- list(...)
    for(. in names(..s)) assign(., ..s[[.]])

    ## labels
    vas <- all.vars(axi)
    ylb <- get0('ilx', inherits=FALSE, ifnotfound=vas[1])
    xlb <- get0('ily', inherits=FALSE, ifnotfound=vas[2])

    ## limites
    vmx <- with(dat, max(val))
    vmi <- with(dat, min(val))
    ylm <- get0('ylim', inherits=FALSE, ifnotfound=c(0, 1.0))

    ## title
    ttl <- get0('itt', inherits=FALSE, ifnotfound=get.unique(dat)$ttl)

    ## box plot
    boxplot(axi, dat, xlab=xlb, ylab=ylb, ylim=ylm, ...)
    abline(1, 0, col='red')
    title(ttl)
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
    if(!is.null(out))
        png(out, width=nc * 3, height=nr * 3, units='in', res=300)

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
