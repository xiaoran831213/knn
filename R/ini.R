## library(MASS)
## library(dplyr, warn.conflicts=FALSE)
## source('R/hlp.R')
## source('R/kp1.R')
## source('R/kd1.R')
## source('R/grm.R')
## source('R/gcta.R')
## source('R/utl.R')
## source('R/lmm.R')
## source("R/mnq/solver.R")


#' @importFrom magrittr %>%
`%>%` <- magrittr::`%>%`

#' @importFrom gtools combinations

#' @useDynLib knn
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function(libpath)
{
    print(paste("Unloading", libpath))
    library.dynam.unload("knn", libpath)
}
