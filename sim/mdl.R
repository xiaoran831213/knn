## models
source('R/utl.R')
source('R/kpl.R')
source('R/ply.R')

ID <- function(x) list(ID=idn(x))
JX <- function(x) list(JX=matrix(1, NROW(x), NROW(x)))
## LN <- function(x) list(LN=ply(scale(x), degree=1))
PL <- function(x) list(PL=ply(x, degree=1))
LP <- function(x) list(LP=lap(x))
## GS <- function(x) list(GS=gau(x))
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
## I2 <- function(x) list(I2=pqw(scale(x), q=2))
## I3 <- function(x) list(I3=pqw(scale(x), q=3))
X2 <- function(x) list(X2=ply(scale(x), degree=2))
X3 <- function(x) list(X3=ply(scale(x), degree=3))


## noise kernel (epsilon)
EPS <- function(x) list(EPS=idn(x))

## othalnormal linear kernels
OL1 <- PK(LN, D=1, O=TRUE)
OL2 <- PK(LN, D=2, O=TRUE)
OL3 <- PK(LN, D=3, O=TRUE)
OL4 <- PK(LN, D=4, O=TRUE)

## linear kernels
## L1 <- PK(LN, D=1, O=FALSE)
## L2 <- PK(LN, D=2, O=FALSE)
## L3 <- PK(LN, D=3, O=FALSE)

## product * Gaussian
OP1 <- PK(PL, D=1, O=TRUE)
OP2 <- PK(PL, D=2, O=TRUE)
OP3 <- PK(PL, D=3, O=TRUE)
PL1 <- PK(PL, D=1, O=FALSE)
PL2 <- PK(PL, D=2, O=FALSE)
PL3 <- PK(PL, D=3, O=FALSE)
GS1 <- PK(GS, D=1, O=FALSE)
GS2 <- PK(GS, D=2, O=FALSE)
GS3 <- PK(GS, D=3, O=FALSE)

## mix
JX1 <- PK(JX, D=1, O=FALSE)
JL1 <- PK(LN, D=1, O=FALSE, J=1)
JL2 <- PK(LN, D=2, O=FALSE, J=2)
JL3 <- PK(LN, D=3, O=FALSE, J=3)
DR1 <- PK(DM, RS, D=1, O=FALSE)
DR2 <- PK(DM, RS, D=2, O=FALSE)
DR3 <- PK(DM, RS, D=3, O=FALSE)

##
JP1 <- JL1
JP2 <- PK(LN, D=2, O=2, J=1)
JP3 <- PK(LN, D=3, O=2, J=1)
JP4 <- PK(LN, D=4, O=2, J=1)
PC2 <- PK(LN, D=2, O=2, J=0)
PC3 <- PK(LN, D=3, O=2, J=0)
PC4 <- PK(LN, D=4, O=2, J=0)


## additive, dominative, recessive
AD1 <- PK(AD, D=1, O=FALSE)
AD2 <- PK(AD, D=2, O=FALSE)
AD3 <- PK(AD, D=3, O=FALSE)
DM1 <- PK(DM, D=1, O=FALSE)
DM2 <- PK(DM, D=2, O=FALSE)
DM3 <- PK(DM, D=3, O=FALSE)
RS1 <- PK(RS, D=1, O=FALSE)
RS2 <- PK(RS, D=2, O=FALSE)
RS3 <- PK(RS, D=3, O=FALSE)
HT1 <- PK(HT, D=1, O=FALSE)
HT2 <- PK(HT, D=2, O=FALSE)
HT3 <- PK(HT, D=3, O=FALSE)
