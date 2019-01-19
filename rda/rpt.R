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
.cp <- function(dat, grp, val='val', cap=0.01)
{
    grp <- split(dat, dat[, grp])
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
        v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
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
    ags <- do.call(rbind, dat %$% 'arg')
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
    rpt <- data.frame(do.call(rbind, dat %$% 'rpt'))
    if(with.arg)
        rpt <- cbind(agg.arg(dvp, ...), rpt)
    rpt
}

agg.par <- function(dvp, cache=TRUE, with.arg=TRUE, drop.fix=FALSE, ...)
{
    dat <- agg.dvp(dvp, cache)
    pas <- dat %$% 'par'
    nms <- unique(unlist(lapply(pas, names)))
    pas <- data.frame(t(sapply(pas, `[`, nms)))
    names(pas) <- nms
    if(drop.fix)
        pas <- pas[, grep('EPS', names(pas), TRUE):ncol(pas)]
    if(with.arg)
        pas <- cbind(agg.arg(dvp, ...), pas)
    pas
}

plt.dvp.par <- function(dvp, out=paste0(dvp, "_dvp_par.png"), cache=TRUE)
{
    ## loading the configurations, estimates, and reports.
    arg <- agg.arg(dvp, cache)
    par <- agg.par(dvp, TRUE, with.arg=FALSE, drop.fix=TRUE)
    rpt <- agg.rpt(dvp, TRUE, with.arg=FALSE)

    ## variance components ralative to sum of squares
    par <- par / rpt$SST

    ## re-arrange
    dat <- cbind(arg, par)
    vas <- c('seed', 'bat', 'mdl', 'mtd', 'pmu', 'std', 'vth') # required
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
    g <- g + facet_grid(mdl~pmu + mtd)

    print(out)
    nfx <- length(unique(dat[, c('pmu', 'mtd', 'vth')]))
    nfy <- length(unique(dat$mdl))
    ufx <- 19 / nfx
    ufy <- 10 / 19 * ufx
    ggsave(out, g, width=19.2, height=min(10, ufy * nfy))
    invisible(dat)
}

plt.dvp.fit <- function(dvp, out=paste0(dvp, "_dvp_fit.png"), cache=TRUE)
{
    ## load reports with configuration, and re-arrange
    dat <- agg.rpt(dvp, cache, with.arg=TRUE)
    vas <- c('seed', 'bat', 'mdl', 'mtd', 'pmu', 'vth') # required
    nxt <- vas[!vas %in% names(dat)]                    # non-exists
    nxt <- matrix(0, nrow(dat), length(nxt), dimnames=list(NULL, nxt))
    dat <- DF(dat, nxt)
    
    library(reshape2)
    dat <- melt(dat, value.name='val', variable.name='key', id.vars=vas)
    
    library(ggplot2)
    dat <- subset(dat, key %in% c('h2b', 'h2c', 'sr2', 'cr2'))
    dat <- .cp(dat, c("pmu", "mtd"), cap=0.05)
    g <- ggplot(dat, aes(x=key, y=val))
    g <- g + geom_boxplot()
    g <- g + coord_cartesian(ylim=c(0, 2.0))
    g <- g + facet_grid(mdl~pmu+vth)

    print(out)
    nfx <- length(unique(dat[, c('pmu', 'mtd')]))
    nfy <- length(unique(dat$mdl))
    ufx <- 19 / nfx
    ufy <- 10 / 19 * ufx
    ggsave(out, g, width=19.2, height=min(10, ufy * nfy))
    invisible(dat)
}

plt.dvp.err <- function(dvp, out=paste0(dvp, "_dvp_err.png"), cache=TRUE)
{
    ## load reports with configuration, and re-arrange
    dat <- agg.rpt(dvp, cache, with.arg=TRUE)
    vas <- c('seed', 'bat', 'mdl', 'mtd', 'pmu', 'std', 'vth') # required
    .nx <- vas[!vas %in% names(dat)]                    # non-exists
    .nx <- matrix(0, nrow(dat), length(.nx), dimnames=list(NULL, .nx))
    dat <- DF(dat, .nx)
    ## dat <- subset(dat, vth < 999)
    dat <- within(dat, vth <- factor(vth))

    library(reshape2)
    dat <- melt(dat, value.name='val', variable.name='key', id.vars=vas)
    dat$pmu <- as.factor(dat$pmu)
    
    library(ggplot2)
    dat <- subset(dat, key %in% c('nlk', 'er2', 'ey2', 'cr1', 'cr2'))
    dat <- .cp(dat, c("pmu", "mtd", "vth"), cap=0.05)
    g <- ggplot(dat, aes(x=vth, y=val))
    ## g <- g + geom_boxplot(aes(color=pmu))
    g <- g + geom_boxplot()
    g <- g + facet_grid(key~mdl, scales='free_y')

    print(out)
    nfx <- length(unique(dat[, c('pmu', 'mtd', 'vth')]))
    nfy <- length(unique(dat$key))    
    ufx <- 19 / nfx
    ufy <- 10 / 19 * ufx
    ggsave(out, g, width=19.2, height=min(10, ufy * nfy))
    invisible(dat)
}

## aggregate reports
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
        arg <- within(arg,
        {
            mdl <- sprintf("%s ~ (%s)^%d", sub("~.*$", "", mdl), sub("~", "", rnd), ply)
            mdl <- gsub(" ", "", mdl)
            rm(pfx, rnd, ply)
        })
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
    agg
}

plt.evl.err <- function(evl, out=paste0(evl, "_evl_err.png"), cap=.05, cache=TRUE)
{
    agg <- agg.evl(evl, cache)
    arg <- with(agg, dat[, arg])
    rpt <- with(agg, dat[, rpt])
    ## rpt <- rpt[, c('sr2', 'sy2', 'nlk', 'ey1', 'er1', 'er2')]
    rpt <- rpt[, c('sr2', 'nlk', 'er1', 'cr1', 'cr2')]
    
    dat <- cbind(arg, rpt)
    vas <- c('mdl', 'nbt', 'bat', 'mtd', 'pmu', 'vth', 'agg') # required
    nxt <- vas[!vas %in% names(dat)]     # non-exists
    nxt <- matrix("def", nrow(dat), length(nxt), dimnames=list(NULL, nxt))
    dat <- DF(dat, nxt)
    
    dat <- subset(dat, cr1 > 0)
    dat <- melt(dat, value.name='val', variable.name='key', id.vars=vas)
    dat <- within(dat, nbt <- as.factor(nbt))
    dat <- within(dat, vth <- as.factor(vth))

    grp <- setdiff(names(dat), c('bat', 'val'))
    if(!is.null(cap))
        dat <- .cp(dat, grp, 'val', cap)
    ## dat <- subset(dat, vth=="Inf" & agg=='mu' & grepl("[\\^]1", mdl))
    dat <- within(dat, tag <- paste(vth, mtd, mdl))
    
    g <- ggplot(dat, aes(x=nbt, y=val))
    g <- g + geom_boxplot(aes(color=agg))
    g <- g + facet_grid(key~tag, scales='free_y')

    print(out)
    nfx <- length(unique(dat[, c('vth', 'mdl')]))
    nfy <- length(unique(dat[, c('key')]))
    ufx <- 19 / nfx
    ufy <- 10 / 19 * ufx
    ggsave(out, g, width=19.2, height=min(10, ufy * nfy))
    rownames(dat) <- NULL
    invisible(dat)
}

plt.evl.mix <- function(evl, out=paste0(evl, "_evl.png"), cap=.05, cache=TRUE)
{
    agg <- agg.evl(evl, cache)
    arg <- with(agg, dat[, arg])
    rpt <- with(agg, dat[, rpt])
    ## kys <- c('ey2', 'h2c')
    kys <- c('ey2', 'cy2')
    dat <- cbind(arg, rpt[, kys])
    dat <- subset(dat, nbt == max(nbt))
    dat <- within(dat,
    {
        nbt <- factor(nbt)
    })
    
    dat <- melt(dat, value.name='val', variable.name='key', measure.vars=kys)
    grp <- setdiff(names(dat), c('bat', 'val'))
    if(!is.null(cap))
    {
        dat <- .cp(dat, grp, 'val', cap)
    }
    dat <- within(dat,
    {
        rsp <- sub("~.*$", "", mdl)
        ply <- sub("^.*[\\^]", "", mdl)
        rm(mdl)
    })
    dat <- subset(dat, rsp %in% c('alc', 'smk'))
    rownames(dat) <- NULL
    
    g <- ggplot(dat, aes(x=rsp, y=val))
    g <- g + geom_boxplot()
    g <- g + facet_wrap(~key, nrow=1L, scales="free", labeller=lbl)
    g <- g + .th
    g <- g + scale_x_discrete(labels=lbl)
        
    print(out)
    nfx <- 1L # length(unique(dat[, c('mdl')]))
    nfy <- length(unique(dat[, c('key')]))
    ## ufx <- 19 / nfx
    ## ufy <- 10 / 19 * ufx
    ## ggsave(out, g, width=19.2, height=ufy * nfy)
    ggsave(out, g, width=19.2, height=10.8/nfy, scale=.6)

    invisible(dat)

}
