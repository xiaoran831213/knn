#' @importFrom magrittr %>%
`%>%` <- magrittr::`%>%`

#' @useDynLib knn
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function(libpath)
{
    print(paste("Unloading", libpath))
    library.dynam.unload("knn", libpath)
}
