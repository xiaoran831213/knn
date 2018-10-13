## methods

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
INM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr=50, cpp=FALSE, psd=0, zbd=0, ...)
## Iterative Z-bound MINQUE
IZM <- function(rsp, knl, ...) knl.mnq(rsp, knl, itr=50, cpp=FALSE, psd=0, zbd=1, ...)

## batched minque
BM0 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^0, ...)
BM1 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^1, ...)
BM2 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^2, ...)
BM3 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^3, ...)
BM4 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^4, ...)
BM5 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^5, ...)
BM6 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^6, ...)
BM7 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^7, ...)
BM8 <- function(rsp, knl, ...) GBT(BMQ, rsp, knl, 2^8, ...)


## GCTA Default
GCT <- function(rsp, knl, ...) gct.rml(rsp, knl, ...)

## MLE
MLE <- function(rsp, knl, ...)
{
    rop.vcm(rsp, knl, ...)
}

## NULL
NUL <- function(rsp, knl, ...)
{
    nul.vcm(rsp, NULL, ...)
}
