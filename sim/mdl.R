## models
source('R/utl.R')
source('R/kpl.R')

ID <- function(x) list(ID=idn(x))
JX <- function(x) list(JX=matrix(1, NROW(x), NROW(x)))
LN <- function(x) list(LN=ply(scale(x), degree=1))
PL <- function(x) list(PL=ply(x, degree=1))
LP <- function(x) list(LP=lap(x))
GS <- function(x) list(GS=gau(x))
IS <- function(x) list(IS=ibs(x))
KN <- function(x) list(KN=kin(x))
RS <- function(x) list(RS=ply((x > 1) * 2 - 1))
DM <- function(x) list(DM=ply((x > 0) * 2 - 1))
AD <- function(x) list(AD=ply((x - 1)))
HT <- function(x) list(HT=ply((x == 1) * 2 - 1))

## genomic models
## A1 <- ~ a
## A2 <- ~ a + I(a^2)
## AA <- ~ a + a:a[]
## AX <- ~ a + I(a^2) + a:a[]
A2 <- function(x) list(A2=p2o(scale(x)))
AA <- function(x) list(AA=p2w(scale(x)))
AX <- function(x) list(AX=ply(scale(x), degree=2))


#' polynomial expansion by kernel
#'
#' ...: the series of basic kernel functions
PK <- function(..., d=2, o=0, j=0)
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
        arg <- c(bas, list(degree=d, coefs=env[['coef']], raw= (o!=1)))
        bas <- do.call(polym, arg)
        
        ## keep orthgonal coeficients
        cfs <- attr(bas, 'coefs')
        env[['coef']] <- if(length(nms) >1) cfs else if(is.null(cfs)) NULL else list(cfs)
        ## print(env[['coef']])

        ## pca?
        if(o == 2)
        {
            pca <- env[['pca']]
            if(is.null(pca))
            {
                print("New PCA.")
                pca <- prcomp(bas, retx=FALSE, scale.=TRUE)
                env[['pca']] <- pca
            }
            else
                print(pca)
            bas <- predict(pca, bas)
            nms <- NULL
        }
        
        bas <- as.data.frame(bas)
        colnames(bas) <- lapply(strsplit(colnames(bas), '[.]'), function(u)
        {
            paste0(nms, u, collapse='.')
        })
        
        bas <- lapply(bas, function(k)
        {
            mtx[utr] <- k; mtx <- t(mtx); mtx[utr] <- k; mtx
        })
        if(j)
        {
            bas <- c(JX(x), bas)
            names(bas)[1] <- 'JX1'
        }
        bas
    }
}

## noise kernel (epsilon)
EPS <- function(x) list(EPS=idn(x))

## othalnormal linear kernels
OL1 <- PK(LN, d=1, o=TRUE)
OL2 <- PK(LN, d=2, o=TRUE)
OL3 <- PK(LN, d=3, o=TRUE)
OL4 <- PK(LN, d=4, o=TRUE)

## linear kernels
LN1 <- PK(LN, d=1, o=FALSE)
LN2 <- PK(LN, d=2, o=FALSE)
LN3 <- PK(LN, d=3, o=FALSE)
LN4 <- PK(LN, d=4, o=FALSE)

## product * Gaussian
OP1 <- PK(PL, d=1, o=TRUE)
OP2 <- PK(PL, d=2, o=TRUE)
OP3 <- PK(PL, d=3, o=TRUE)
PL1 <- PK(PL, d=1, o=FALSE)
PL2 <- PK(PL, d=2, o=FALSE)
PL3 <- PK(PL, d=3, o=FALSE)
GS1 <- PK(GS, d=1, o=FALSE)
GS2 <- PK(GS, d=2, o=FALSE)
GS3 <- PK(GS, d=3, o=FALSE)

## mix
JX1 <- PK(JX, d=1, o=FALSE)
JL1 <- PK(LN, d=1, o=FALSE, j=1)
JL2 <- PK(LN, d=2, o=FALSE, j=2)
JL3 <- PK(LN, d=3, o=FALSE, j=3)
DR1 <- PK(DM, RS, d=1, o=FALSE)
DR2 <- PK(DM, RS, d=2, o=FALSE)
DR3 <- PK(DM, RS, d=3, o=FALSE)

##
JP1 <- JL1
JP2 <- PK(LN, d=2, o=2, j=1)
JP3 <- PK(LN, d=3, o=2, j=1)
JP4 <- PK(LN, d=4, o=2, j=1)
PC2 <- PK(LN, d=2, o=2, j=0)
PC3 <- PK(LN, d=3, o=2, j=0)
PC4 <- PK(LN, d=4, o=2, j=0)


## additive, dominative, recessive
AD1 <- PK(AD, d=1, o=FALSE)
AD2 <- PK(AD, d=2, o=FALSE)
AD3 <- PK(AD, d=3, o=FALSE)
DM1 <- PK(DM, d=1, o=FALSE)
DM2 <- PK(DM, d=2, o=FALSE)
DM3 <- PK(DM, d=3, o=FALSE)
RS1 <- PK(RS, d=1, o=FALSE)
RS2 <- PK(RS, d=2, o=FALSE)
RS3 <- PK(RS, d=3, o=FALSE)
HT1 <- PK(HT, d=1, o=FALSE)
HT2 <- PK(HT, d=2, o=FALSE)
HT3 <- PK(HT, d=3, o=FALSE)

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
