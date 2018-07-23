## Training Track
msg <- function(obj, fmt=NULL, hdr=NULL, ...)
{
    if(is.function(obj))
    {
        ret <- function()
        {
            sprintf(fmt, obj(...))
        }
    }
    else if(inherits(obj, 'formula'))
    {
        ret <- function()
        {
            e <- environment(obj)
            sprintf(fmt, eval(attr(terms(obj), 'variables'), e)[[1]])
        }
    }
    else
    {
        ret <- function()
        {
            sprintf(fmt, obj)
        }
    }

    if(is.null(fmt))
        fmt <- "%8s"

    if(is.null(hdr))
        hdr <- deparse(substitute(obj))
    hdr <- sub("^~", "", hdr)
    
    len <- as.integer(sub("^%[+-]?([0-9]*)[.]?([0-9])?[A-z]$", "\\1", fmt))
    hdr <- sprintf(paste0("% ", len, "s"), hdr)
    structure(ret, class=c('msg', 'function'), hdr=hdr)
}

## is
is.msg <- function(.) 'msg' %in% class(.)

## Header of the tracks
hd.msg <- function(...)
{
    d <- Filter(is.msg, unlist(list(...)))
    d <- sapply(d, function(.)
    {
        h <- attr(., "hdr")
    })
    paste(d, collapse=" ")
}

## A line of the tracks
ln.msg <- function(...)
{
    d <- Filter(is.msg, unlist(list(...)))
    d <- sapply(d, do.call, list())
    paste(d, collapse=" ")
}
