source('R/hlp.R')                       # helpers
source('R/utl.R')                       # utilities
library(dplyr, warn.conflicts=FALSE)
library(ggplot2, warn.conflicts=FALSE)
library(reshape2)

## shared theme
.th <- theme(
    axis.title.x=element_blank(), axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    strip.text.x = element_text(size=12, face="bold"),
    strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="red", fill="#CCCCFF"),
    legend.title=element_blank(), legend.position='bottom')

## cap the values
.cp <- function(dat, grp, val='val', cap=0.01, mtd=c('both', 'upper', 'lower'))
{
    grp <- split(dat, dat[, grp])
    mtd <- match.arg(mtd, c('both', 'upper', 'lower'))
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        if(mtd == 'upper')
            v <- pmin(v, quantile(v, 1-cap, na.rm=TRUE))
        else if(mtd == 'lower')
            v <- pmax(v, quantile(v, 0+cap, na.rm=TRUE))
        else
        {
            v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
            v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
        }
        g[, val] <- v
        g
    })
    dat <- do.call(rbind, grp)
    dat
}

## labeler
lbl <- function (labels, multi_line = TRUE)
{
    labels <- label_value(labels, multi_line = multi_line)
    dc <- c(
        alc="Alcohol", smk="Smoking", cnb="Cannabis",
        mse="MSE", cyh="r", rsq="r^2", hsq="h^2", nlk="NLK",
        sr2="rSLP[2]",
        er2="MSE(g, hat(g))",
        cr2="COR(g, hat(g))",
        ey2="MSE(y, hat(y))",
        cy2="COR(y, hat(y))",
        h2c="h^2")
    if (multi_line)
    {
        lapply(unname(labels), lapply, function(values)
        {
            if(values %in% names(dc))
                values <- dc[[values]]
            c(parse(text = as.character(values)))
        })
    }
    else
    {
        lapply(labels, function(values)
        {
            values <- paste0("list(", values, ")")
            lapply(values, function(expr) c(parse(text = expr)))
        })
    }
}

## aggregate model development
agg.dvp <- function(dvp, cache=TRUE)
{
    rds <- paste0(dvp, '.rds')
    if(cache && file.exists(rds))
    {
        rps <- readRDS(rds)
    }
    else
    {
        rps <- dir(dvp, '[.]rpt$', full=TRUE)
        rps <- lapply(rps, function(f)
        {
            print(f)
            readRDS(f)
        })
        saveRDS(rps, rds)
    }
    invisible(rps)
}


agg.arg <- function(dvp, cache=TRUE, simple=TRUE)
{
    dat <- agg.dvp(dvp, cache)
    ags <- do.call(.rbd, dat %$% 'arg')
    if(is.null(ags$ply))
        ags$ply <- 1
    if(simple)
    {
        ags <- within(ags,
        {
            mdl <- sprintf("%s ~ (%s)^%d", sub("~.*$", "", mdl), sub("~", "", rnd), ply)
            mdl <- gsub(" ", "", mdl)
            rm(pfx, rnd, ply)
        })
    }
    ags
}


agg.rpt <- function(dvp, cache=TRUE, with.arg=TRUE, ...)
{
    dat <- agg.dvp(dvp, cache)
    rpt <- data.frame(do.call(.rbd, dat %$% 'rpt'))
    if(with.arg)
        rpt <- cbind(agg.arg(dvp, ...), rpt)
    rpt
}


agg.par <- function(dvp, cache=TRUE, with.arg=TRUE, drop.fix=FALSE, ...)
{
    dat <- agg.dvp(dvp, cache)
    ## pas <- dat %$% 'par'
    ## nms <- unique(unlist(lapply(pas, names)))
    ## pas <- data.frame(t(sapply(pas, `[`, nms)))
    ## names(pas) <- nms
    pas <- do.call(.rbd, dat %$% 'par')
    if(drop.fix)
        pas <- pas[, grep('EPS', names(pas), TRUE):ncol(pas)]
    if(with.arg)
        pas <- cbind(agg.arg(dvp, ...), pas)
    pas
}


plt.par <- function(dvp, out=paste0(dvp, "_par.png"), cache=TRUE)
{
    ## loading the configurations, estimates, and reports.
    arg <- agg.arg(dvp, cache)
    par <- agg.par(dvp, TRUE, with.arg=FALSE, drop.fix=TRUE)
    rpt <- agg.rpt(dvp, TRUE, with.arg=FALSE)

    ## variance components ralative to sum of squares
    par <- par / rpt$SST

    ## re-arrange
    dat <- cbind(arg, par)
    vas <- c('seed', 'bat', 'mdl', 'mtd', 'pmu', 'tag', 'msl') # required
    nxt <- vas[!vas %in% names(dat)]                    # non-exists
    nxt <- matrix(0, nrow(dat), length(nxt), dimnames=list(NULL, nxt))
    dat <- DF(dat, nxt)
    
    library(reshape2)
    dat <- melt(dat, value.name='val', variable.name='key', id.vars=vas)
    dat <- subset(dat, !is.na(val))
    rownames(dat) <- NULL
    
    library(ggplot2)
    g <- ggplot(dat, aes(x=key, y=val))
    g <- g + geom_boxplot()
    g <- g + coord_cartesian(ylim=c(-.2, 1.2))
    g <- g + facet_grid(tag~msl)

    print(out)
    nfy <- length(unique(dat$tag))
    ufy <- 10 / nfy
    nfx <- length(unique(dat$msl))
    ufx <- 19 / nfx
    if(ufx / ufy < 19 / 10)
        ufy <- ufx / 19 * 10
    else
        ufx <- ufy / 10 * 19
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy))
    invisible(dat)
}

plt.fit <- function(dvp, out=paste0(dvp, "_fit.png"), cache=TRUE)
{
    ## load reports with configuration, and re-arrange
    dat <- agg.rpt(dvp, cache, with.arg=TRUE)
    vas <- c('seed', 'bat', 'mdl', 'mtd', 'pmu', 'tag', 'msl') # required
    nxt <- vas[!vas %in% names(dat)]                    # non-exists
    nxt <- matrix(0, nrow(dat), length(nxt), dimnames=list(NULL, nxt))
    dat <- DF(dat, nxt)
    
    library(reshape2)
    dat <- melt(dat, value.name='val', variable.name='key', id.vars=vas)
    
    library(ggplot2)
    dat <- subset(dat, key %in% c('hsq', 'yel', 'ycl', 'zel', 'zcl', 'nlk'))
    dat <- .cp(dat, c("mtd", "tag", "key"), cap=0.08, mtd='u')
    g <- ggplot(dat, aes(x=tag, y=val))
    g <- g + geom_boxplot()
    ## g <- g + coord_cartesian(ylim=c(0, 2.0))
    g <- g + facet_grid(key~msl, scales="free")

    print(out)
    nfy <- length(unique(dat$key))
    ufy <- 10 / nfy
    nfx <- length(unique(dat$msl))
    ufx <- 19 / nfx
    if(ufx / ufy < 19 / 10)
        ufy <- ufx / 19 * 10
    else
        ufx <- ufy / 10 * 19
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy))
    invisible(dat)
}


## aggregate model eveluation
agg.evl <- function(evl, cache=TRUE)
{
    rds <- paste0(evl, '.rds')
    if(cache && file.exists(rds))
        agg <- readRDS(rds)
    else
    {
        agg <- lapply(dir(evl, '[.]rds$', full=TRUE), function(f)
        {
            print(f)
            r <- readRDS(f)
            if(inherits(r, 'try-error'))
                r <- NULL
            r
        })
        rpt <- do.call(.rbd, agg %$% 'rpt')
        arg <- do.call(.rbd, agg %$% 'arg')
        arg$mdl <- NULL
        arg$rnd <- NULL
        par <- do.call(.rbd, agg %$% 'par')
        dat <- cbind(arg, par, rpt)
        
        rpt <- colnames(rpt)
        arg <- colnames(arg)
        par <- colnames(par)

        vcs <- par[grep('EPS', par, TRUE):length(par)]
        fix <- par[!par %in% vcs]
        agg <- list(dat=dat, arg=arg, par=par, rpt=rpt, vcs=vcs, fix=fix)
        saveRDS(agg, rds)
    }
    invisible(agg)
}

plt.evl <- function(evl, out=paste0(evl, "_evl.png"), cap=.05, cache=TRUE)
{
    agg <- agg.evl(evl, cache)
    arg <- with(agg, dat[, arg])
    rpt <- with(agg, dat[, rpt])
    kys <- c('hsq', 'yel', 'ycl', 'nlk')
    rpt <- rpt[, kys]
    dat <- DF(cbind(arg, rpt))

    dat <- melt(dat, value.name='val', variable.name='key', measure.vars=kys)
    dat <- within(dat, nbt <- as.factor(nbt))
    ## dat <- within(dat, msl <- as.factor(msl))
    
    ## grp <- setdiff(names(dat), c('bat', 'val', 'msl'))
    ## if(!is.null(cap))
    ##     dat <- .cp(dat, grp, 'val', cap)
    ## dat <- subset(dat, agg=='mu')
    
    g <- ggplot(dat, aes(x=tag, y=val))
    ## g <- g + geom_boxplot(aes(color=agg))
    g <- g + geom_boxplot()
    g <- g + facet_grid(key~mtd, scales='free_y')

    print(out)
    nfy <- length(unique(dat$key))
    ufy <- 10 / nfy
    nfx <- length(unique(dat$mtd))
    ufx <- 19 / nfx
    if(ufx / ufy < 19 / 10)
        ufy <- ufx / 19 * 10
    else
        ufx <- ufy / 10 * 19
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy))
    rownames(dat) <- NULL
    invisible(dat)
}

plt.ref <- function(evl, out=paste0(evl, "_evl.png"), cap=.05, cache=TRUE)
{
    agg <- agg.evl(evl, cache)
    arg <- with(agg, dat[, arg])
    rpt <- with(agg, dat[, rpt])
    ## kys <- c('hsq', 'yel', 'ycl', 'nlk')
    kys <- c('yel', 'ycl')
    rpt <- rpt[, kys]
    dat <- DF(cbind(arg, rpt))
    ## dat <- subset(dat, tag != 'fw2c')
    
    grp <- dat[, c('bat', 'pfx', 'evl')]
    ref <- by(dat, grp, function(g)
    {
        g <- within(g,
        {
            yel <- yel - yel[tag == 'wgs1']
            ## nlk <- nlk - nlk[tag == 'wgs1']
            ycl <- ycl - ycl[tag == 'wgs1']
        })
        g
    })
    dat <- do.call(rbind, ref)

    dat <- melt(dat, value.name='val', variable.name='key', measure.vars=kys)
    dat <- .cp(dat, 'tag', 'val', cap,  "lower")
 
    g <- ggplot(dat, aes(x=tag, y=val))
    g <- g + geom_boxplot(aes(color=evl))
    g <- g + facet_grid(key~pfx, scales='free_y')

    print(out)
    nfy <- length(unique(dat$key))
    ufy <- 10 / nfy
    nfx <- length(unique(dat$pfx))
    ufx <- 19 / nfx
    if(ufx / ufy < 19 / 10)
        ufy <- ufx / 19 * 10
    else
        ufx <- ufy / 10 * 19
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy))
    rownames(dat) <- NULL
    invisible(dat)
}

main.rpt <- function()
{
    agg.dvp('rda/dvp/alc_ck2_nlk', FALSE)
    agg.dvp('rda/dvp/alc_ck3_nlk', FALSE)
    agg.dvp('rda/dvp/sbp_ck2_nlk', FALSE)
    agg.dvp('rda/dvp/sbp_ck3_nlk', FALSE)
    agg.dvp('rda/dvp/smk_ck2_nlk', FALSE)
    agg.dvp('rda/dvp/smk_ck3_nlk', FALSE)

    agg.evl('rda/evl/alc_ck2_nlk', FALSE)
    agg.evl('rda/evl/alc_ck3_nlk', FALSE)
    agg.evl('rda/evl/sbp_ck2_nlk', FALSE)
    agg.evl('rda/evl/sbp_ck3_nlk', FALSE)
    agg.evl('rda/evl/smk_ck2_nlk', FALSE)
    agg.evl('rda/evl/smk_ck3_nlk', FALSE)

    plt.ref('rda/evl/alc_ck2_nlk')
    plt.ref('rda/evl/alc_ck3_nlk')
    plt.ref('rda/evl/sbp_ck2_nlk')
    plt.ref('rda/evl/sbp_ck3_nlk')
    plt.ref('rda/evl/smk_ck2_nlk')
    plt.ref('rda/evl/smk_ck3_nlk')
}
