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


#' regular expression substitution
#'
#' This simplifies the syntax of regular expression matching.
#'
#' @param text text to match against pattern
#' @param pat regular expression, and replacement
#' @return a vector of the sam elength as text, with substituted items
#' @examples
#' > text <- c('a','b','c')
#' > text %s% c('[ab]','_')
#' [1] "_" "_" "c"
`%s%` <- function(text, pat){
	if (length(pat)>1){
		rep <- pat[2]
	} else {
		rep <- ""
	}
	return(gsub(pattern=pat[1],replacement=rep,text))
}


#' printf prints
#'
#' a wrapper around sprintf but it actually outputs on screen.
#' @param ... arguments to sprintf
#' @return nothing
printf <- function(...){
 cat(sprintf(...),sep='')
}
