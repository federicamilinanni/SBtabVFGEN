#' Make C compatible names
#'
#' Uses make.names internally, but replaces dots with underscores.
#'
#' @param Labels a character vector with words that need to be turned
#'	 into names.
#' @export
#' @return a vector with unique names, can be used as variable names
#'	 in C
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
#' @export
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
#' [1] "a"	 "b"	 "c"	 "d"	 "1"	 "1 / 2"
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
#'	 v
#' @return a matrix with various versions of v (columns) one per
#'	 setting described in data.frame Table. The names can have a ">"
#'	 prefix in the names (see SBtab rules)
#' @export
#' @examples
#' > v<-c(1,2,3)
#' > names(v)<-c('a','b','c')
#'
#' > data<-data.frame(row.names=c('low','med','high'),
#'				  b=c(0.5,2.5,5.5),
#'				  comment=c('b < 1','close to default','b > 2×default'))
#'
#'		b		  comment
#' low  0.5			b < 1
#' med  2.5 close to default
#' high 5.5	b > 2×default
#'
#' > update_from_table(v,data)
#'   low med high
#' a   1   1	1
#' b   2   2	2
#' c   3   3	3
update_from_table <- function(v,Table){
	if (is.null(v) || is.null(Table)) return(NULL)
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
#'	 (either)
#' @return a vector with names corresponding to !ID and values taken
#'	 from !DefaultValue, !Value, or !InitialValue
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

#' A parser for .ods files with SBtab document structure
#'
#' This uses the readODS package to read the file.
#' SBtab sheets are themselves tables dedicated to a specific type of
#' model property: Reaction, Compound, Parameter, etc.
#'
#' The SBtab content is not interpreted in any way.
#' @param ods.file a string (file's name)
#' @return SBtab
#' a list of tables (data.frames), one per ods sheet
#' SBtab[['TableName']] retrives a data.frame
#' comment(SBtab) is the name of the document
#'
#' @keywords import
#' @examples
#' model.sbtab<-sbtab_from_ods('model.ods')
#' @export
sbtab_from_ods <- function(ods.file){
	M <- readODS::read.ods(ods.file)
	lM <- length(M)
	SBtab <- list(length=lM)
	header <- paste0(M[[1]][1,],collapse='')
	document.name <- sbtab.header.value(header,'Document')
	table.name <- vector(length=lM)
	for (i in 1:lM){
		header <- paste0(M[[i]][1,],collapse='')
		table.name[i] <- sbtab.header.value(header,'TableName')
		SBtab[[i]] <- M[[i]][-c(1,2),]
		names(SBtab[[i]]) <- M[[i]][2,]
	}
	names(SBtab) <- table.name
	cat("Tables found: ",paste(table.name,collapse=', '),"\n")
	comment(SBtab) <- document.name
	return(SBtab)
}

#' A parser for a bunch of .tsv files with SBtab document content
#'
#' This function reads the files (not using any package).
#' SBtab sheets are themselves tables dedicated to a specific type of
#' model property: Reaction, Compound, Parameter, etc.
#'
#' The SBtab content is not interpreted in any way.
#' @param tsv.file a character vector (file names, one per sheet),
#'	 defaults to all tsv files in the current directory.
#' @return SBtab a list of tables, one per file in tsv.file list
#'	 SBtab[['Reaction']] retrieves the table of reactions, a
#'	 data.frame comment(SBtab) retrieves the SBtab document name
#' @keywords import
#' @examples
#' model.sbtab<-sbtab_from_tsv(dir(pattern='.*[.]tsv$'))
#' @export
sbtab_from_tsv <- function(tsv.file=dir(pattern='[.]tsv$')){
	SBtab <- list()
	header <- readLines(tsv.file[1],n=1)
	document.name <- sbtab.header.value(header,"Document")
	printf("[tsv] file[1] «%s» belongs to Document «%s»\n\tI'll take this as the Model Name.\n",tsv.file[1],document.name)
	for (f in tsv.file){
		header <- readLines(f,n=1)
		TableName <- sbtab.header.value(header,'TableName')
		SBtab[[TableName]] <- read.delim(f,as.is=TRUE,skip=1,check.names=FALSE,comment.char="%",blank.lines.skip=TRUE)
		attr(SBtab[[TableName]],"TableName") <- TableName
	}
	comment(SBtab) <- document.name
	return(SBtab)
}

#' collect all variable names
#'
#' this function collects all varibale names in SBtab, here we use the
#' !ID column (because the !Name column can be something that isn't a
#' programming language kind of variable-name).
#'
#' If the react flag is FALSE then really all names are collected,
#' otherwise only those that can appear in reaction formulae and
#' kinetic laws of reactions.
#'
#' @param tab a list of lists with the table content
#' @param react if TRUE, this function only collects the IDs of
#'	 Compouds, Parameters and Expressions
#' @return a vector of all names
all.vars <- function(tab,react=TRUE){
	allvars <- c(tab$Parameter[["!ID"]],tab$Compound[["!ID"]])
	tNames <- names(tab)
	if (react){
		l <- tNames %in% c("Input","Constant","Expression")
		tNames <- tNames[l]
	}
	for (T in tNames){
		if ("!ID" %in% names(tab[[T]])) {
			allvars <- c(allvars,tab[[T]][["!ID"]])
		}
	}
	return(trimws(unlist(allvars)))
}

#' Checks that all SBtab references are valid
#'
#' A reference is a column name that starts with ">": >Calcium is a
#' reference to a compound (probably) called Calcium (or a Parameter,
#' etc.)
#'
#' @param tab a list of tables
#' @param allvars a list of all variable IDs
all.refs.valid <- function(tab,allvars=all.vars(tab,reac=FALSE)){
	r <- TRUE
	for (T in tab){
		cNames <- colnames(T)
		l <- grepl("^[~>][a-zA-Z]\\w*$",cNames)
		refs <- cNames[l] %s% "[~>]"
		valid.refs <- refs %in% allvars
		if (!all(valid.refs)){
			printf("Table %s contains invalid references:\n",attr(T,"TableName"))
			printf("%30s\n",refs[!valid.refs])
			r <- FALSE
			warning("Not all refs are valid.")
		}
	}
	return(r)
}

#' Are the IDs the same as the Names
#'
#' This is not mandatory, but there could be errors if somewhere in
#' the document a Name is used inside of mathematical expressions
#' while the ID was needed. This can help find those mistakes.
#'
#' @param T one table (data.frame)
id.eq.name <- function(T){
	if (!all(c("!ID","!Name") %in% names(T))) return(NA)
	ID <- T[["!ID"]]
	Name <- T[["!Name"]]
	l <- ID == Name
	if (any(!l)){
		cat("these IDs and Names are different:\n")
		printf("%30s    %s\n","ID","Name")
		printf("%30s != %s\n",ID[!l],Name[!l])
	}
	return(l)
}

#' validate spelling in SBtab document
#'
#' This function checks the SBtab document for consistency.
#' It finds misspelled varibale names in the reaction kinetics.
#'
#' @param tab a list of lists as returned by `sbtab_from_tsv()`
sbtab.mistakes <- function(tab){
	stopifnot("Reaction" %in% names(tab))
	vars <- ftsplit(gsub("[-()+*/]+"," ",tab$Reaction[["!KineticLaw"]]),"[ ]+",re=TRUE)
	av <- all.vars(tab)
	l <- vars %in% av
	if (any(!l)) {
		cat("These variables appear in the reaction table, but are not defined in the rest of the document\n")
		print(vars[!l])
	}
	## point out where ID and Name are different, maybe that is a problem
	for (T in tab){
		l <- id.eq.name(T)
		if (!is.na(l) && !all(l)){
			warning("Not all IDs are equal to the Name attribute, but maybe this is on purpose.")
		}
	}
	if (all.refs.valid(tab)){
		cat("All internal references are valid (~REF and >REF)\n")
	}
}

#' Read data from SBtab
#'
#' This function reads all datasets from SBtab Documents that are
#' referenced in the Experiment table and returns them as a list of
#' data.frames
#'
#' The returned list always describes simulation experiments:
#' time series experiments, possibly with scheduled events.
sbtab.data <- function(tab){
	l <- grepl("([tT]able([oO]f)?)?[eE]xperiments?|[mM]easurements?|[sS]imulations?",names(tab))
	if (any(l) && sum(l)==1) {
		E <- tab %1% l
	} else {
		warning("There must be exactly one Table of Experiments.")
		return(NA)
	}
	n <- dim(E)[1]
	experiments <- vector("list",length=n)
	names(experiments) <- E[["!ID"]]
	l <- grepl("!([eE]xperiment)?Type",names(E))
	if (any(l)){
		type <- E %1% l
		time.series <- grepl("time[_ ]series",type,ignore.case=TRUE)
		dose.response <- grepl("dose[_ ]response",type,ignore.case=TRUE)
		dose.sequence <- grepl("dose[ _]sequence",type,ignore.case=TRUE)
	} else {
		time.series <- rep(TRUE,n) # by default
		dose.response <- !time.series
		dose.sequence <- !time.series
	}
	v <- tab$Compound[["!InitialValue"]]
	names(v) <- tab$Compound[["!ID"]]
	initVal <- update_from_table(v,E)

	v <- tab$Input[["!DefaultValue"]]
	names(v) <- tab$Input[["!ID"]]
	input <- update_from_table(v,E)

	id <- E[["!ID"]]
	for (i in 1:n){
		stopifnot(id[i] %in% names(tab))
		tNames <- names(tab[[id[i]]])
		if (time.series[i]){
			l <- grepl(">([a-zA-Z][^ ]*)",tNames)
			outNames <- tNames[l]
			errNames <- tNames[l] %s% c(">","~")
			experiments[[id[i]]][["outputValues"]] <- tab[[id[i]]][l]
			experiments[[id[i]]][["errorValues"]] <- tab[[id[i]]][errNames]
			names(experiments[[id[i]]][["outputValues"]]) <- outNames %s% ">"
			names(experiments[[id[i]]][["errorValues"]]) <- errNames %s% "~"
			experiments[[id[i]]][["input"]] <- input[,i]
			experiments[[id[i]]][["initialState"]] <- initVal[,i]
			experiments[[id[i]]][["outputTimes"]] <- tab[[id[i]]][["!Time"]]
		}
	}
	return(experiments)
}
