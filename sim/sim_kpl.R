source('R/kpl.R')                       # kernel players

source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")
source("R/utl.R")
source("sim/utl.R")
kpl.t1 <- function(N=100, P=300)
{
    x <- readRDS('data/1kg_c05.rds')[1:N, 1:P]

    k1 <- kin(x)
    k2 <- gau(x)
    k3 <- esn(x, p=1)
    k4 <- esn(x, p=.5)
    
    par(mfrow=c(2, 2))
    image(k1)
    image(k2)
    image(k3)
    image(k4)
}
