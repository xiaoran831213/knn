library(dplyr)

readRpt <- function(fdr, ..., ref=FALSE)
{
    rpt <- list()
    for(f in dir(c(fdr, ...), '.rds$', full.names=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next

        ## use null model as a measuring stick
        r <- subset(r, dat == 'evl', -dat)
        mse <- subset(r, key == 'mse' & mtd != 'nul')
        nul <- subset(r, key == 'mse' & mtd == 'nul', val)
        mse <- within(mse, val <- val / unlist(nul))

        nlk <- subset(r, key == 'nlk' & mtd != 'nul')
        nul <- subset(r, key == 'nlk' & mtd == 'nul', val)
        nlk <- within(nlk, val <- val / unlist(nul))
        
        r <- rbind(mse, nlk)
        rpt <- c(rpt, list(r))
    }
    rpt <- do.call(rbind, rpt)
    rpt <- within(rpt,
    {
        mat <- sub('[.].*$', '', mtd)
    })

    clp <- function(x)
    {
        x <- gsub("^c", "", x)
        x <- gsub("[ ()]", "", x)
        x <- gsub(",", "+", x)
        x
    }
    rpt <- within(rpt,
    {
        sim <- sprintf("%s(%s)~%s", lnk, clp(oks), clp(yks))
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
    ## vq1 <- with(dat, quantile(val, .95))
    vq0 <- with(dat, quantile(val, .05))
    ylm <- get0('ylim', inherits=FALSE, ifnotfound=c(0, 1.0))

    ## title
    ttl <- get.unique(dat)$ttl

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
            ipt(gp[[gc]], ...)
            gc <- gc + 1
        }
    }
    title(ttl, outer=TRUE, adj=0)

    if(!is.null(out))
        dev.off()
}

m1 <- function()
{
    s0 <- subset(readRpt('sim/run/s00', ref=1), key=='nlk', -c(key))
    s1 <- subset(readRpt('sim/run/s01', ref=1), key=='nlk', -c(key))
    s2 <- subset(readRpt('sim/run/s02', ref=1), key=='nlk', -c(key))
    s3 <- subset(readRpt('sim/run/s03', ref=1), key=='nlk', -c(key))

    plt(s0, oxy=.~het, ipt=bxp)
    plt(s1, oxy=.~het, ipt=bxp)
    plt(s2, oxy=.~het, ipt=bxp)
    plt(s3, oxy=.~het, ipt=bxp)

    da <- rbind(s0, s1, s2, s3)
    d1 <- within(subset(da, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    d2 <- within(subset(da, grepl('[.]rop$', mtd)), {use <- 'rop'; mtd <- sub('[.]rop$', '', mtd)})
    
    plt(d1, oxy=sim~het, ipt=bxp, axi=val~mtd, out='rpt/2018_07_08/km1_mnq.png')
    plt(d2, oxy=sim~het, ipt=bxp, axi=val~mtd, out='rpt/2018_07_08/km1_rop.png')

    plt(da, oxy=sim~het, ipt=bxp, axi=val~mtd)
}

m2 <- function()
{
    s0 <- subset(readRpt('sim/run/s10', ref=1), key=='nlk', -c(key))
    s1 <- subset(readRpt('sim/run/s11', ref=1), key=='nlk', -c(key))
    s2 <- subset(readRpt('sim/run/s12', ref=1), key=='nlk', -c(key))
    s3 <- subset(readRpt('sim/run/s13', ref=1), key=='nlk', -c(key))
    da <- rbind(s0, s1, s2, s3)

    plt(s0, oxy=.~het, ipt=bxp)
    plt(s1, oxy=.~het, ipt=bxp)
    plt(s2, oxy=.~het, ipt=bxp)
    plt(s3, oxy=.~het, ipt=bxp)

    d1 <- within(subset(da, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    d2 <- within(subset(da, grepl('[.]rop$', mtd)), {use <- 'rop'; mtd <- sub('[.]rop$', '', mtd)})
    
    plt(d1, oxy=sim~het, ipt=bxp, axi=val~mtd, out='rpt/2018_07_08/km2_mnq.png')
    plt(d2, oxy=sim~het, ipt=bxp, axi=val~mtd, out='rpt/2018_07_08/km2_rop.png')
}
