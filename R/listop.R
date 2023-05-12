#' first valid item in data object
#'
#' Returns the first item in list/data.frame x that is TRUE according
#' to l, where l is a logical vector. This is similar to dplyr::first
#' The vector l describes some previously established condition for
#' the list.
#'
#' if x is a data.frame, return the first column that
#' is labeled as TRUE by l.
#'
#' @param x list or data.frame
#' @param l the result of some condition (appropriate for x)
#' @return the first item in x that has a TRUE l value
#' @examples
#' > x <- list(a=data.frame(x=c(0,1),y=c(100,128)),b=data.frame(x=c(1,2),y=c(128,200)))
#' > l <- grepl('[bB]',names(x))
#' > m <- grepl('[yY]',names(x[[1]]))
#' > x %1% l %1% m
#' [1] 128 200
`%1%`  <- function(x,l){
	stopifnot(is.logical(l))
	if (is.list(x[l])) {
		return(x[l][[1]])
	} else if(is.data.frame(x[l])){
		return(x[l][,1])
	}
}
