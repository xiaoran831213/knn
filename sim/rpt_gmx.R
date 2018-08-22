source('sim/rpt.R')

m0 <- function(name='b20')
{
    fdr <- file.path('sim/run', name)
    err <- readRpt(fdr)
    rtm <- readRtm(fdr)
    dat <- within(rbind(rtm, err), {H <- N * R; N <- N * Q})
    dat <- subset(dat, se=-c(Q, R))

    img <- file.path('rpt/2018_08_23/img', name)
    unlink(img, TRUE, TRUE)
    for(m in c('mnq', 'mle'))           # method
    {
        for(a in c('cyh', 'loo', 'mse', 'nlk', 'ssz')) # aggregation
        {
            . <- c('gct', paste0(m, '.', c('avg', 'whl', 'ssz', a)))
            for(k in c('cyh', 'mse', 'nlk', 'loo'))
            {
                dt <- subset(dat, mtd %in% . & key == k)
                dt <- within(dt, mtd <- sub('^.*[.]', '', mtd))
                od <- file.path(img, k)

                dir.create(od, FALSE, TRUE)
                fn <- file.path(od, paste0('mbt_', m, '_', a, '.png'))
                print(fn)
                plt(dt, oxy=bsz~sim, ipt=bxp, out=fn)
            }
        }

        ## running time
        {
            . <- c('gct', paste0(m, '.', c('bat', 'whl')))
            k <- 'rtm'
            {
                dt <- subset(dat, mtd %in% . & key == k)
                dt <- within(dt, mtd <- sub('^.*[.]', '', mtd))
                dt <- within(dt, mtd <- relevel(as.factor(mtd), 'gct'))
                od <- file.path(img, k)

                dir.create(od, FALSE, TRUE)
                fn <- file.path(od, paste0('mbt_', m, '.png'))
                print(fn)
                plt(dt, oxy=bsz~sim, ipt=bxp, out=fn)
            }
        }
    }
}

nlk <- function(name='c10')
{
    library(ggplot2)
    f <- file.path('sim/run', name)
    d <- readRpt(f)
    d <- subset(d, key=='nlk', -c(key))
    d <- within(d, mtd <- sub('[.].*$', '', mtd))

    for(p in unique(d$P))
    {
        g <- ggplot(subset(d, P==p), aes(x=mtd, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(eps ~ sc)
        f <- paste0(file.path('~', name), '_', p, '.png')
        ggsave(f, g, width=10, height=10)
    }
}
