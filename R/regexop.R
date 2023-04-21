#' regular sub-expression matching
#'
#' This simplifies the syntax of regular expression matching.
#'
#' @param text text to match against pattern
#' @param pattern regular expression pattern (can have subexpressions)
#' @return a list of the same length as text, with lists of matches in each item
#' @examples match some math expressions:
#' x <- c('2*x^2','x','3')
#' x %~% '([0-9]*)[*]?x?(\\^[0-9]+)?'
#' [[1]]
#' [1] "2*x^2" "2"     "^2"
#'
#' [[2]]
#' [1] "x" ""  ""
#'
#' [[3]]
#' [1] "3" "3" ""
`%~%` <- function(text, pattern){
	r <- regexec(pattern=pattern,text=text)
	m <- regmatches(text,r)
	return(m)
}
