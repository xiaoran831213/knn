## methods

## batched minque
BM1 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 100, ...)
BM2 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 200, ...)
BM3 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 300, ...)
BM4 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 400, ...)
BM5 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 500, ...)
BM6 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 600, ...)
BM7 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 700, ...)
BM8 <- function(rsp, knl, ...) GBT(knl.mnq, rsp, knl, 800, ...)

## PSD minque
PMQ <- function(rsp, knl, ...)
{
    knl.mnq(rsp, knl, psd=TRUE, ...)
}

## Normal MINQUe
NMQ <- function(rsp, knl, ...)
{
    knl.mnq(rsp, knl, psd=FALSE, ...)
}

## GCTA Default
GCT <- function(rsp, knl, ...)
{
    gct.rml(rsp, knl, ...)
}

## MLE
MLE <- function(rsp, knl, ...)
{
    rop.vcm(rsp, knl, ...)
}

## NULL
NUL <- function(rsp, knl, ...)
{
    nul.vcm(rsp, knl, ...)
}
