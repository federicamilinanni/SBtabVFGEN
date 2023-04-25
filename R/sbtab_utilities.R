#' Make C compatible names
#'
#' Uses make.names internally, but replaces dots with underscores.
#'
#' @param Labels a character vector with words that need to be turned
#'     into names.
#' @return a vector with unique names, can be used as variable names
#'     in C
make.cnames <- function(Labels){
	Names <- gsub("'([^']*)'","lsquo\\1rsquo",trimws(Labels))
	Names <- gsub('"([^"]*)"',"ldquo\\1rdquo",Names)
	Names <- gsub("&","and",Names,fixed=TRUE)
	Names <- gsub("|","or",Names,fixed=TRUE)
	Names <- gsub("'","prime",Names)
	Unique.Names <- gsub("[.]","_",make.names(Names, unique = TRUE, allow_ = TRUE))
	return(Unique.Names)
}

#' Fixed-token-split a scalar string into parts
#'
#' This function splits a string like strsplit, but it removes
#' whitespace from the result (on both ends) and the split token is
#' taken literally
#'
#' @param str a string (character vector of length 1)
#' @param s a split token
#' @param re defaults to FALSE, if TRUE s is treated as a regular expression
#' @return a character vector of the components without leading or trailing whitespace
#' @examples ftsplit(" A + 2*B ","+")
#' [1] "A"   "2*B"
#' @examples x<-c('a+b','c+d','1 + 2'); lapply(x,ftsplit,'+')
#' [[1]]
#' [1] "a" "b"
#' [[2]]
#' [1] "c" "d"
#' [[3]]
#' [1] "1" "2"
#' @examples this also works, but mixes up the components:
#' x<-c('a+b','c+d+1','1 / 2'); ftsplit(x,'+')
#' [1] "a"     "b"     "c"     "d"     "1"     "1 / 2"
ftsplit <- function(str,s,re=FALSE){
	return(trimws(unlist(strsplit(str,s,fixed=!re))))
}

#' Update a named vector with values from a table
#'
#' This function can be used in the circumstance that you already have
#' a default vector and want to update some of its entries from a
#' data.frame with more specific values. The table describes the
#' circumstamces of these more specific values and has columns named
#' like the vector elements. The function returns a matrix with one
#' column per scenario, updated using the table. This is used while
#' parsing SBtab tables.
#'
#' The returned value M has several copies of vector v
#' with some values changed accoring to a table (Qantity Matrix) If
#' the table has columns that name members of v (also named) they will
#' be used. The M will have as many columns as the table has rows.
#'
#' @param v the vector to update
#' @param Table a table with column names partially matching those in
#'     v
#' @return a matrix with various versions of v (columns) one per
#'     setting described in data.frame Table. The names can have a ">"
#'     prefix in the names (see SBtab rules)
#' @export
#' @examples
#' > v<-c(1,2,3)
#' > names(v)<-c('a','b','c')
#'
#' > data<-data.frame(row.names=c('low','med','high'),
#'                  b=c(0.5,2.5,5.5),
#'                  comment=c('b < 1','close to default','b > 2×default'))
#'
#'        b          comment
#' low  0.5            b < 1
#' med  2.5 close to default
#' high 5.5    b > 2×default
#'
#' > update_from_table(v,data)
#'   low med high
#' a   1   1    1
#' b   2   2    2
#' c   3   3    3
update_from_table <- function(v,Table){
    N <- names(v)
    stopifnot(is.data.frame(Table))
    n <- nrow(Table)
    M <- matrix(v,nrow=length(v),ncol=n)
    rownames(M)<-names(v)
    l <- startsWith(names(Table),'>')
    T<-Table[l]
    names(T)<-sub(">","",names(Table[l]))
    NT <- N[N %in% names(T)]
    M[NT,] <- t(T[NT])
    colnames(M)<-rownames(Table)
    return(M)
}

#' read vector from sbtab Quantity table
#'
#' This function extracts a vector from an sbtab table, and names the
#' elements according to the "!ID" column.
#'
#' @param Table the sbtab table imported via sbtab_from_{ods,tsv}
#'     (either)
#' @return a vector with names corresponding to !ID and values taken
#'     from !DefaultValue, !Value, or !InitialValue
#' @keywords sbtab quantity
#' @export
sbtab_quantity <- function(Table){
    N <- names(Table)
    NT <- N[N %in% c('!DefaultValue','!Value','!InitialValue')]
    v <- unlist(Table[NT])
    n <- unlist(Table['!ID'])
    names(v) <- n
    return(v)
}

#' Read a property from the SBtab header
#'
#' This function retrieves a property from an sbtab header:
#' PropertyName='PropertyValue' (properties are these key=value
#' pairs).
#'
#' @param sbtab.header the first line of an SBtab file
#' @param key the left side of each key=value pair
#' @return value of the key=value pair
sbtab.header.value <- function(sbtab.header,key='Document'){
	## the table title has to be in either one of the columns of row 1
	m <- unlist(sbtab.header %~% sprintf("%s='([^']+)'",key))
	if (length(m)>0){
		property <- m[2] # so the first experssion in parentheses
	} else {
		stop(sprintf("property «%s» not set in SBtab header: «%s».",key,header))
	}
	return(property)
}

#' printf prints
#'
#' a wrapper around sprintf but it actually outputs on screen.
#' @param ... arguments to sprintf
#' @return nothing
printf <- function(...){
 cat(sprintf(...),sep='')
}
