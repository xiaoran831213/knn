source('sim/rpt.R')

m1 <- function()
{
    rds <- 'sim/run/b1x.rds'
    if(file.exists(rds))
    {
        d <- readRDS(rds)
    }
    else
    {
        b0 <- readRpt('sim/run/b10')
        b1 <- readRpt('sim/run/b11')
        b2 <- readRpt('sim/run/b12')
        b3 <- readRpt('sim/run/b13')
        d <- rbind(b0, b1, b2, b3)
        saveRDS(d, rds)
    }
    d <- within(d, {H <- N * R; N <- N * Q})
    d <- subset(d, se=-c(Q, R))

    f1 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        x <- subset(d, sim=='I(ga)~p2')
        x <- subset(x, mtd %in% c('mle.cyh', 'mle.whl', 'mnq.cyh', 'mnq.whl', 'gct'))
        plt(subset(x, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_nlk_cyh.png')
        plt(subset(x, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_mse_cyh.png')
    }

    f2 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        x <- subset(d, sim=='I(ga)~p2')
        x <- subset(x, mtd %in% c('mle.ssz', 'mle.whl', 'mnq.ssz', 'mnq.whl', 'gct'))
        plt(subset(x, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_nlk_ssz.png')
        plt(subset(x, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_mse_ssz.png')
    }

    g1 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(p2)~p2')
        d <- subset(d, mtd %in% c('mle.cyh', 'mle.whl', 'mnq.cyh', 'mnq.whl', 'gct'))
        plt(subset(d, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_nlk_cyh.png')
        plt(subset(d, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_mse_cyh.png')
    }

    g2 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(p2)~p2')
        d <- subset(d, mtd %in% c('mle.ssz', 'mle.whl', 'mnq.ssz', 'mnq.whl', 'gct'))
        plt(subset(d, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_nlk_ssz.png')
        plt(subset(d, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_mse_ssz.png')
    }

    f3 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(ga)~p2')
        d <- subset(d, mtd %in% c('mle.bat', 'mle.whl', 'mnq.bat', 'mnq.whl', 'gct'))
        plt(subset(d, key=='rtm'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_rtm.png')
    }
}

m4 <- function()
{
    rds <- 'sim/run/b4x.rds'
    if(file.exists(rds))
    {
        d <- readRDS(rds)
    }
    else
    {
        b0 <- readRpt('sim/run/b40')
        b1 <- readRpt('sim/run/b41')
        b2 <- readRpt('sim/run/b42')
        b3 <- readRpt('sim/run/b43')
        d <- rbind(b0, b1, b2, b3)
        saveRDS(d, rds)
    }
    d <- within(d, {H <- N * R; N <- N * Q})
    d <- subset(d, se=-c(Q, R))

    f1 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(ga)~p1')
        d <- subset(d, mtd %in% c('mle.cyh', 'mle.whl', 'mnq.cyh', 'mnq.whl', 'gct'))
        plt(subset(d, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_nlk_cyh.png')
        plt(subset(d, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_mse_cyh.png')
    }

    f2 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(ga)~p1')
        d <- subset(d, mtd %in% c('mle.ssz', 'mle.whl', 'mnq.ssz', 'mnq.whl', 'gct'))
        plt(subset(d, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_nlk_ssz.png')
        plt(subset(d, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_mse_ssz.png')
    }

    g1 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(p2)~p2')
        d <- subset(d, mtd %in% c('mle.cyh', 'mle.whl', 'mnq.cyh', 'mnq.whl', 'gct'))
        plt(subset(d, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_nlk_cyh.png')
        plt(subset(d, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_mse_cyh.png')
    }

    g2 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(p2)~p2')
        d <- subset(d, mtd %in% c('mle.ssz', 'mle.whl', 'mnq.ssz', 'mnq.whl', 'gct'))
        plt(subset(d, key=='nlk'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_nlk_ssz.png')
        plt(subset(d, key=='mse'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_mse_ssz.png')
    }

    f3 <- function()
    {
        ## kmq mle.cyh mle.loo mle.mse mle.nlk mle.ssz mle.whl     mnq mnq.ssz mnq.whl 
        d <- subset(d, sim=='I(ga)~p2')
        d <- subset(d, mtd %in% c('mle.bat', 'mle.whl', 'mnq.bat', 'mnq.whl', 'gct'))
        plt(subset(d, key=='rtm'), oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt4_rtm.png')
    }
}

c0 <- function()
{
    d <- readRpt('sim/run/c00', 'sim/run/c01', 'sim/run/c02', 'sim/run/c03')
}
