source('sim/rpt.R')

m0 <- function()
{
    b0 <- readRpt('sim/run/b00', ref=1)
    b1 <- readRpt('sim/run/b01', ref=1)
    b2 <- readRpt('sim/run/b02', ref=1)

    b0.mse <- subset(b0, key=='mse', -key)
    plt(b0.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt0_mse.png')
    b0.nlk <- subset(b0, key=='nlk', -key)
    plt(b0.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt0_nlk.png')
    
    b1.mse <- subset(b1, key=='mse', -key)
    plt(b1.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_mse.png')
    b1.nlk <- subset(b1, key=='nlk', -key)
    plt(b1.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt1_nlk.png')

    b2.mse <- subset(b2, key=='mse', -key)
    plt(b2.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt2_mse.png')
    b2.nlk <- subset(b2, key=='nlk', -key)
    plt(b2.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/bt2_nlk.png')
}

m1 <- function()
{
    b0 <- readRpt('sim/run/b10', ref=1)
    b1 <- readRpt('sim/run/b11', ref=1)
    b2 <- readRpt('sim/run/b12', ref=1)

    b0.mse <- subset(b0, key=='mse', -key)
    plt(b0.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b10_mse.png')
    b0.nlk <- subset(b0, key=='nlk', -key)
    plt(b0.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b10_nlk.png')
    
    b1.mse <- subset(b1, key=='mse', -key)
    plt(b1.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b11_mse.png')
    b1.nlk <- subset(b1, key=='nlk', -key)
    plt(b1.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b11_nlk.png')

    b2.mse <- subset(b2, key=='mse', -key)
    plt(b2.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b12_mse.png')
    b2.nlk <- subset(b2, key=='nlk', -key)
    plt(b2.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b12_nlk.png')
}

m2 <- function()
{
    b0 <- readRpt('sim/run/b20')
    b1 <- readRpt('sim/run/b21')
    b2 <- readRpt('sim/run/b22')

    b0 <- subset(b0, het %in% c(0.0, 0.5))
    b1 <- subset(b1, het %in% c(0.0, 0.5))
    b2 <- subset(b2, het %in% c(0.0, 0.5))

    b0 <- subset(b0, mtd != 'rop')
    b1 <- subset(b1, mtd != 'rop')
    b2 <- subset(b2, mtd != 'rop')

    b0.mse <- subset(b0, key=='mse', -key)
    plt(b0.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b20_mse.png')
    b0.nlk <- subset(b0, key=='nlk', -key)
    plt(b0.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b20_nlk.png')
    
    b1.mse <- subset(b1, key=='mse', -key)
    plt(b1.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b21_mse.png')
    b1.nlk <- subset(b1, key=='nlk', -key)
    plt(b1.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b21_nlk.png')

    b2.mse <- subset(b2, key=='mse', -key)
    plt(b2.mse, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b22_mse.png')
    b2.nlk <- subset(b2, key=='nlk', -key)
    plt(b2.nlk, oxy=het~bsz, ipt=bxp, out='rpt/2018_07_15/b22_nlk.png')
}

m3 <- function()
{
    rds <- 'sim/run/b3x.rds'
    if(file.exists(rds))
        da <- readRDS(rds)
    else
    {
        b0 <- readRpt('sim/run/b30', ref=1)
        b1 <- readRpt('sim/run/b31', ref=1)
        b2 <- readRpt('sim/run/b32', ref=1)
        b3 <- readRpt('sim/run/b33', ref=1)
        da <- rbind(b0, b1, b2, b3)
        saveRDS(da, rds)
    }
    da <- within(da, {H <- N * R; N <- N * Q})
    da <- subset(da, se=-c(Q, R))
    
    ## take out MNQUE and GCTA
    d1 <- subset(da, sim=="st(ga)~p1" & mtd=="gct", -sim) # gcta
    d2 <- subset(da, sim=="st(ga)~p1" & mtd=="mnq", -sim) # minque 1
    d3 <- subset(da, sim=="st(ga)~p2" & mtd=="mnq", -sim) # minque 2
    d2 <- within(d2, mtd <- "mq1")
    d3 <- within(d3, mtd <- "mq2")
    db <- rbind(d1, d2, d3)

    dc <- subset(db, het==0.0, -het)
    plt(dc, oxy=.~key, ipt=bxp, axi=val~mtd, out='rpt/kma/img/jsm_vcm_bmk.png')
}
