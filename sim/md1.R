## models
source('R/utl.R')
source('R/kpl.R')

ID <- function(x) list(ID=diag(1, NROW(x)))
JX <- function(x) list(JX=matrix(1, NROW(x), NROW(x)))
LN <- function(x) list(LN=ply(scale(x), degree=1))
PL <- function(x) list(PL=ply(x, degree=1))
LP <- function(x) list(LP=lap(x))
GS <- function(x) list(GS=gau(x))
IS <- function(x) list(IS=ibs(x))
KN <- function(x) list(KN=kin(x))
RS <- function(x)
{
    x <- (x > 1) * 2 - 1
    x <- x[, apply(x, 2L, sd) > 0]
    x <- scale(x)
    list(RS=ply(x))
}
DM <- function(x)
{
    x <- (x > 0) * 2 - 1
    x <- x[, apply(x, 2L, sd) > 0]
    x <- scale(x)
    list(DM=ply(x))
}
AD <- function(x)
{
    x <- x - 1
    x <- scale(x)
    list(AD=ply(x))
}
HT <- function(x)
{
    x <- (x == 1) * 2 - 1
    x <- x[, apply(x, 2L, sd) > 0]
    x <- scale(x)
    list(HT=ply(x))
}

#' polynomial expansion by kernel
#'
#' ...: the series of basic kernel functions
PK <- function(..., d=2, orth=FALSE, J=FALSE)
{
    env <- environment()
    function(x)
    {
        ## print(env[['coef']])
        N <- NROW(x)                    # sample size
        mtx <- matrix(0, N, N)          # sample matrix
        utr <- upper.tri(mtx, TRUE)     # sample upper.tri

        bas <- c(...)                   # bases
        bas <- unlist(lapply(bas, do.call, list(x=x)), FALSE)
        bas <- lapply(bas, `[`, utr)
        nms <- names(bas)
        
        ## expansion
        arg <- c(bas, list(degree=d, coefs=env[['coef']], raw=!orth))
        bas <- do.call(polym, arg)
        
        ## keep orthgonal coeficients
        cfs <- attr(bas, 'coefs')
        env[['coef']] <- if(length(nms) >1) cfs else if(is.null(cfs)) NULL else list(cfs)
        ## print(env[['coef']])

        bas <- as.data.frame(bas)
        colnames(bas) <- lapply(strsplit(colnames(bas), '[.]'), function(u)
        {
            paste0(nms, u, collapse='.')
        })
        
        lapply(bas, function(k)
        {
            mtx[utr] <- k; mtx <- t(mtx); mtx[utr] <- k; mtx
        })
    }
}

## othalnormal linear kernels
OL2 <- PK(LN, d=2, orth=TRUE)
OL3 <- PK(LN, d=3, orth=TRUE)

## linear kernels
LN1 <- LN
LN2 <- PK(LN, d=2, orth=FALSE)
LN3 <- PK(LN, d=3, orth=FALSE)
LN4 <- PK(LN, d=4, orth=FALSE)

## product * Gaussian
OP2 <- PK(PL, d=2, orth=TRUE)
OP3 <- PK(PL, d=3, orth=TRUE)
PL1 <- PL
PL2 <- PK(PL, d=2, orth=FALSE)
PL3 <- PK(PL, d=3, orth=FALSE)

GS1 <- GS

## mix
JMX <- JX
DR1 <- PK(DM, RS, d=1, orth=FALSE)
DR2 <- PK(DM, RS, d=2, orth=FALSE)
DR3 <- PK(DM, RS, d=3, orth=FALSE)
DR4 <- PK(DM, RS, d=4, orth=FALSE)

## additive, dominative, recessive
AD1 <- PK(AD, d=1, orth=FALSE)
AD2 <- PK(AD, d=2, orth=FALSE)
AD3 <- PK(AD, d=3, orth=FALSE)
DM1 <- PK(DM, d=1, orth=FALSE)
DM2 <- PK(DM, d=2, orth=FALSE)
DM3 <- PK(DM, d=3, orth=FALSE)
RS1 <- PK(RS, d=1, orth=FALSE)
RS2 <- PK(RS, d=2, orth=FALSE)
RS3 <- PK(RS, d=3, orth=FALSE)
HT1 <- PK(HT, d=1, orth=FALSE)
HT2 <- PK(HT, d=2, orth=FALSE)
HT3 <- PK(HT, d=3, orth=FALSE)

coef <- function(k)
{
    env <- environment(x)
    if(!is.null(env))
        environment(k)[['coef']]
    else
        NULL
}
`coef<-` <- function(x, value)
{
    env <- environment(x)
    if(!is.null(env))
        env[['coef']] <- value
    x
}
