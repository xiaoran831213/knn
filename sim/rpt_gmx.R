source('sim/rpt.R')

m0 <- function()
{
    rds <- 'sim/run/b0x.rds'
    if(file.exists(rds))
    {
        d <- readRDS(rds)
    }
    else
    {
        . <- paste0('sim/run/b0', seq(0, 3))
        err <- readRpt(.)
        rtm <- readRtm(.)
        d <- rbind(rtm, err)
        d <- within(d, {H <- N * R; N <- N * Q})
        d <- subset(d, se=-c(Q, R))
        saveRDS(d, rds)
    }

    d1 <- subset(d, sim=='i2(p2)~p2' & key != 'rtm')
    ## d <- subset(d, het %in% c(0.0, 0.5, 1.0))

    for(m in c('mnq', 'mle'))           # method
    {
        for(a in c('cyh', 'loo', 'mse', 'nlk')) # aggregation
        {
            . <- c('gct', paste0(m, '.', c('avg', 'whl', 'ssz', a)))
            print(.)

            ## dev.new()
            dt <- within(subset(d1, mtd %in% .), mtd <- sub('^.*[.]', '', mtd))
            fn <- paste0('rpt/2018_08_09/img/mbt_', m, '_', a, '_mse.png')
            plt(subset(dt, key=='mse'), oxy=bsz~het, ipt=bxp, out=fn)
            fn <- paste0('rpt/2018_08_09/img/mbt_', m, '_', a, '_nlk.png')
            plt(subset(dt, key=='nlk'), oxy=bsz~het, ipt=bxp, out=fn)
        }
    }
}

c0 <- function()
{
    d <- readRpt('sim/run/c00', 'sim/run/c01', 'sim/run/c02', 'sim/run/c03')
}
