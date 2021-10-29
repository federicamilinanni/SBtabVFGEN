#' copy vector n times and change the columns according to (n×m) Table
#' 
#' Update a vector v to a matrix M with several copies of that vector
#' with some values changed accoring to a table (Qantity Matrix) If
#' the table has columns that name members of v (also named) they will
#' be used. The M will have as many columns as the table has rows.
#'
#' @param v a vector with named elements
#' @param Table an sbtab table with columns that reference members of
#'     v via the ">" prefix (e.g. the column ">ATP" will change the
#'     v['ATP'] value)
#' @return M, a matrix with n columns, where n is the number of rows
#'     in Table.
#' @keywords sbtab quantity matrix
#' @export
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
