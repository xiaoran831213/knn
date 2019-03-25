## methods
library(wrapGCTA)
library(mnq)

## batched MLE
BL0 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^0, ...)
BL1 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^1, ...)
BL2 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^2, ...)
BL3 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^3, ...)
BL4 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^4, ...)
BL5 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^5, ...)
BL6 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^6, ...)
BL7 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^7, ...)
BL8 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^8, ...)
BL9 <- function(rsp, knl, ...) GBT(rop.vcm, rsp, knl, 2^9, ...)

## Ordinary PSD      MINQUE
OPM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr= 1, cpp=FALSE, psd=1, zbd=0, ...)
## ordinary Normal   MINQU
ONM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr= 1, cpp=FALSE, psd=0, zbd=0, ...)
## ordinary Z-bound  MINQUE
OZM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr= 1, cpp=FALSE, psd=0, zbd=1, ...)
## Iterative PSD     MINQUE
IPM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr=50, cpp=FALSE, psd=1, zbd=0, ...)
## Iterative Normal  MINQUE
INM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr=50, cpp=FALSE, psd=0, ebd=1, vbd=1, zbd=0, ...)
## Iterative Z-bound MINQUE
IZM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr=50, cpp=FALSE, psd=0, ebd=1, vbd=1, zbd=1, ...)

MNQ <- function(rsp, kns, xmx=NULL, ...) fwd(rsp, kns, xmx, tol=1e-4, rpt=1, ...)

## batched minque
UBZ <- 32
BM0 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^0, ...)
BM1 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^1, ...)
BM2 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^2, ...)
BM3 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^3, ...)
BM4 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^4, ...)
BM5 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^5, ...)
BM6 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^6, ...)
BM7 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^7, ...)
BM8 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^8, ...)
BM9 <- function(rsp, knl, xmx=NULL, ...) GBT(MNQ, rsp, knl, xmx, UBZ * 2^9, ...)


MQ0 <- function(y, K, itr=20, ...) mnq(y, K, zbd=0, vbd=0, ebd=0, itr=itr, ...)
MQ1 <- function(y, K, itr=20, ...) mnq(y, K, zbd=0, vbd=1, ebd=0, itr=itr, ...)
MQ2 <- function(y, K, itr=20, ...) mnq(y, K, zbd=0, vbd=0, ebd=1, itr=itr, ...)
MQ3 <- function(y, K, itr=20, ...) mnq(y, K, zbd=0, vbd=1, ebd=1, itr=itr, ...)

ZQ0 <- function(y, K, itr=20, ...) knl.mnq(y, K, zbd=1, vbd=0, ebd=0, itr=itr, ...)
ZQ1 <- function(y, K, itr=20, ...) knl.mnq(y, K, zbd=1, vbd=1, ebd=0, itr=itr, ...)
ZQ2 <- function(y, K, itr=20, ...) knl.mnq(y, K, zbd=1, vbd=0, ebd=1, itr=itr, ...)
ZQ3 <- function(y, K, itr=20, ...) knl.mnq(y, K, zbd=1, vbd=1, ebd=1, itr=itr, ...)

ZGC <- function(y, K, itr=100, ...) gcta.reml(y, K, itr=itr, zbd=1, alg=1)
UGC <- function(y, K, itr=100, ...) gcta.reml(y, K, itr=itr, zbd=0, alg=1)

## GCTA Default
GCT <- function(rsp, knl, ...)
{
    r <- gcta.reml(rsp, knl, quiet=FALSE, ...)
    par <- r$par
    rpt <- vpd(rsp, knl, w=r$par, ...)
    list(par=par, rpt=rpt, rtm=c(rtm=r$rtm))
}


## MLE
MLE <- function(rsp, knl, ...)
{
    rop.vcm(rsp, knl, ...)
}

## NULL
NUL <- function(rsp, knl, ...)
{
    nul.vcm(rsp, ...)
}
