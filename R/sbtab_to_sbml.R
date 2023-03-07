.unit.kind <- function(kind){
	stopifnot(is.character(kind) && length(kind)==1)
	if (grepl("^(l|L|litre|liter)$",kind)){
		k <- "litre"
	} else if (grepl("^(mole?)$",kind)) {
		k <- "mole"
	} else if (grepl("^(m|meter|metre)$",kind)) {
		k <- "metre"
	} else if (grepl("^(kg|kilogram)$",kind)){
		k <- "kilogram"
	} else if (grepl("^(g|gram)$",kind)){
		k <- "gram"
	} else if (grepl("^(A|ampere)$",kind)){
		k <- "ampere"
	} else if (grepl("^(cd|candela)$",kind)){
		k <- "candela"
	} else if (grepl("^(s|second)$",kind)){
		k <- "second"
	} else if (grepl("^(K|kelvin)$",kind)){
		k <- "kelvin"
	} else {
		k <- "dimensionless"
	}
	return(k)
}

.unit.scale <- function(prefix){
	stopifnot(is.character(prefix) && length(prefix)==1)
	if (grepl("^G$|^giga$",prefix)){
		s <- 9
	} else if (grepl("^M$|^mega$",prefix)){
		s <- 6
	} else if (grepl("^k$|^kilo$",prefix)){
		s <- 3
	} else if (grepl("^h$|^hecto$",prefix)){
		s <- 2
	} else if (grepl("^d$|^deci$",prefix)){
		s <- -1
	} else if (grepl("^c$|^centi$",prefix)){
		s <- -2
	} else if (grepl("^m$|^milli$",prefix)){
		s <- -3
	} else if (grepl("^u$|^µ$|^micro$",prefix)){
		s <- -6
	} else if (grepl("^n$|^nano$",prefix)){
		s <- -9
	} else if (grepl("^p$|^pico$",prefix)){
		s <- -12
	} else if (grepl("^f$|^femto$",prefix)){
		s <- -15
	} else {
		s <- 0
	}
	return(s)
}

.unit.id.from.string <- function(unit.str,prnt=FALSE){
	uid <- unit.str
	uid <- sub("^1$","dimensionless",uid)
	uid <- gsub("1/","one_over_",uid)
	uid <- gsub("/","_per_",uid)
	uid <- gsub("[*[:blank:]]","_",uid)
	uid <- gsub("[()]","",uid)
	uid <- gsub("\\^2","_square",uid)
	uid <- gsub("\\^3","_cube",uid)
	uid <- gsub("\\^([0-9]+)","_to_the_power_of_\1",uid)
	uid <- make.names(uid,unique=FALSE)
	if (prnt){
		message("units in «!Unit» column:")
		print(unit.str)
		message("automatically created sbml unit ids:")
		print(uid)
	}
	return(uid)
}

#' Unit Interpreter
#'
#' This function will try its best to interpret strings like
#'	"liter/(nmol ms)"
#' rules: 1. only one slash is allowed
#'        2. M can be mega or mol/l
#'        3. prefixes and units can be words or single letters
#'        4. everything after a slash is the denominator
#'        5. u is an accepted replacement for μ
#'        6. no parentheses
#' this retruns a data.frame with components as in the sbml standard:
#' kind, multiplier, scale and exponent
#' since there is only one slash,parentheses do nothing
#' everything after a slash is the denominator, so: l/mol s is the same as (l)/(mol s)
#' Remark: not all units are understood.
#' @param unit.str a string that contains a human readable unit
#' @return data.frame with an interpretation of the unit
.interpret.unit.from.string <- function(unit.str){
	stopifnot(length(unit.str)==1)
	.kind <- NULL
	.multiplier <- NULL
	.scale <- NULL
	.exponent <- NULL
	message(sprintf("---%11s---",unit.str))
	print(unit.str)
	a <- gsub("[()]","",unit.str)
	a <- gsub("molarity","mol l^-1",a);
	if (grepl("/",unit.str)){
		a <- unlist(strsplit(a,split="/",fixed=TRUE))
		message(sprintf("«%s» is interpreted as:\n\tNumerator «%s»\n\tDenominator: «%s»\n",unit.str,a[1],a[2]))
	}
	n <- length(a)
	stopifnot(n==1 || n==2)

	for (j in 1:n){
		b <- unlist(strsplit(a[j],split="[* ]"))
		for (u in b){
			pat <- paste0("^(G|giga|M|mega|k|kilo|c|centi|m|milli|u|μ|micro|n|nano|p|pico|f|femto)?",
			              "(l|L|liter|litre|g|gram|mole?|s|second|m|meter|metre|K|kelvin|cd|candela|A|ampere)",
			              "\\^?([-+]?[0-9]+)?$")
			##print(pat)
			r <- regexec(pattern=pat,text=u) #,perl=TRUE)

			if (u == "1"){
			    .u.s <- 0
			    .u.k <- "dimensionless"
			    x <- 1
			} else if (r[[1]][1] > 0){
			    m <- unlist(regmatches(u, r))
			    print(m)
			    .u.s <- .unit.scale(m[2])
			    .u.k <- .unit.kind(m[3])
			    if (nchar(m[4])>0){
			        x <- switch(j,as.numeric(m[4]),-as.numeric(m[4]))
			    }else{
			        x <- switch(j,1,-1)
			    }
			} else {
			    warning(sprintf("unit «%s» didn't match any pre-defined unit, because of «%s». (According to the limited [dumb] rules defined in this R script)",unit.str,u))
			    .u.s <- 0
			    .u.k <- "dimensionless"
			    x <- 1
			}
			.multiplier <- c(.multiplier,1.0)
			.scale <- c(.scale,.u.s)
			.kind <- c(.kind,.u.k)
			.exponent <- c(.exponent,x)

			message(sprintf("«%s» has been interpreted as: (%s×10^(%i))^(%i)",u,.u.k,.u.s,x))
		}
	}
	unit <- data.frame(scale=.scale,multiplier=.multiplier,exponent=.exponent,kind=.kind)
	print(unit)
	message(sprintf("---%-11s---","done"))

	return(unit)
}


.create.sbml.unit <- function(sbml,unit,id){
	n <- nrow(unit)
	unitdef <- Model_createUnitDefinition(sbml);
	UnitDefinition_setId(unitdef,id);
	for (i in 1:n){
		u <- UnitDefinition_createUnit(unitdef);
		Kind <- switch(unit$kind[i],litre="UNIT_KIND_LITRE",metre="UNIT_KIND_METRE",second="UNIT_KIND_SECOND",mole="UNIT_KIND_MOLE",gram="UNIT_KIND_GRAM","UNIT_KIND_DIMENSIONLESS")
		Unit_setKind(u, Kind);
		Unit_setExponent(u,unit$exponent[i]);
		Unit_setMultiplier(u,unit$multiplier[i]);
		Unit_setScale(u,unit$scale[i]);
	}
}

.unique.units <- function(sbml,Unit){
	U <- unique(Unit)
	uid <- .unit.id.from.string(U,prnt=TRUE)
	unit <- list()
	for (i in 1:length(uid)){
		unit[[i]] <- .interpret.unit.from.string(U[i])
		.create.sbml.unit(sbml,unit[[i]],uid[i])
	}
	names(unit) <- uid
	return(unit)
}


.sbtab.parameter.to.sbml <- function(sbml,Parameter){
	all.uid <- .unit.id.from.string(Parameter$Unit)
	num.parameters <- nrow(Parameter)
	Name <- row.names(Parameter)

	if ("Value" %in% names(Parameter)){
		Value <- as.numeric(Parameter$Value)
	} else if ("DefaultValue" %in% names(Parameter)) {
		Value <- as.numeric(Parameter$DefaultValue)
	} else {
		stop("Parameters have no Value column")
	}
	for (i in 1:num.parameters){
		message(sprintf("adding SBtab parameter %i: «%s»",i,Name[i]))
		para <- Model_createParameter(sbml)
		Parameter_setId(para,Name[i])
		Parameter_setName(para,Name[i])

		this.unit.id <- all.uid[i]
		message(sprintf("unit of parameter %i (%s): «%s»",i,Name[i],this.unit.id))
		Parameter_setUnits(para,this.unit.id)
		Parameter_setValue(para,Value[i])
	}
}

.sbtab.compound.to.sbml <- function(sbml,Compound,CompName,SubstanceUnitID){
	all.uid <- .unit.id.from.string(Compound$Unit)
	num.species <- nrow(Compound)
	Name <- row.names(Compound)
	for (i in 1:num.species){
		print(Name[i])
		sp <- Model_createSpecies(sbml)
		this.unit.id <- all.uid[i]
		message(sprintf("unit of species %i (%s): «%s»",i,Name[i],this.unit.id))
		Species_setId(sp, Name[i])

		Species_setCompartment(sp, CompName)
		Species_setUnits(sp,SubstanceUnitID)
		Species_setInitialConcentration(sp, Compound$InitialValue[i])
		Ai <- Compound$Assignment[i]
		if (Compound$IsConstant[i]){
			Species_setBoundaryCondition(sp,"true")
			Species_setName(sp, Name[i])
		} else if (!grepl("^(0|FALSE|NONE|NULL)?$",Ai)){
			Species_setBoundaryCondition(sp,"true")
			Species_setName(sp, Ai)
		} else {
			Species_setName(sp, Name[i])
		}
	}
}

.sbtab.constant.to.sbml <- function(sbml,Constant){
	all.uid <- .unit.id.from.string(Constant$Unit)
	num.constants <- nrow(Constant)
	print(num.constants)
	Name <- row.names(Constant)
	for (i in 1:num.constants){
		message(sprintf("adding SBtab Constant %i: «%s»",i,Name[i]))
		print(Name[i])
		const <- Model_createParameter(sbml)
		Parameter_setId(const,Name[i])
		Parameter_setName(const,Name[i])
		this.unit.id <- all.uid[i]
		message(sprintf("unit of constant %i (%s): «%s»",i,Name[i],this.unit.id))
		Parameter_setUnits(const, this.unit.id)
		Parameter_setValue(const, as.numeric(Constant$Value[i]))
	}
}

.sbtab.reaction.to.sbml <- function(sbml,Reaction,Compartment){
	num.reactions <- nrow(Reaction)
	Name <- row.names(Reaction)
	OneCompartmentModel=FALSE;
	if (is.na(Compartment) || length(Compartment$ID==1)){
		OneCompartmentModel=TRUE;
	}
	AB <- strsplit(Reaction$Formula,split="<=>")
	for (i in 1:num.reactions){
		message(sprintf("adding reaction %i: «%s»",i,Name[i]))
		reaction <- Model_createReaction(sbml);

		Reaction_setId(reaction,Name[i]);
		Reaction_setName(reaction,Name[i]);
		Reaction_setReversible(reaction, Reaction$IsReversible[i]);

		message(sprintf("converting flux to mathml: «%s»",Reaction$Flux[i]))
		kl <- Reaction_createKineticLaw(reaction);
		Flux <- Reaction$Flux[i]
		if (OneCompartmentModel){
			Flux <- sprintf("(%s)*%s",Flux,Compartment$ID[1]);
			message("reaction rates in SBML need to have the unit «substance/time» rather than «concentration/time».\nSo conventional kinetics need to be multiplied by compartment volumes.");
			message("Since this model has only one Compartment, all kinetics will be multiplied by its volume for the SBML model.");
		}
		astMath <- parseL3Formula(Flux); # or: readMathMLFromString(mathXMLString);
		KineticLaw_setMath(kl, astMath);

		A <- trimws(strsplit(AB[[i]][1],"+",fixed=TRUE)[[1]])
		l <- grepl("^\\s*(NULL|NIL|Ø)\\s*$",A)
		A <- A[!l]
		B <- trimws(strsplit(AB[[i]][2],"+",fixed=TRUE)[[1]])
		l <- grepl("^\\s*(NULL|NIL|Ø)\\s*$",B)
		B <- B[!l]

		message("Reactants: ")
		print(A)
		message("Products: ")
		print(B)
		message("---")
		for (a in A){
			r<-regexec(pattern="([0-9]*)\\s*[*]?\\s*(\\w+)",text=a)
			m <- unlist(regmatches(a,r))
			print(m)
			RefName <- m[3]
			spr  <-  Reaction_createReactant(reaction)
			SimpleSpeciesReference_setSpecies(spr, RefName)
			if (nchar(m[2])>0){
			    SpeciesReference_setStoichiometry(spr,as.numeric(m[2]))
			}
		}
		for (b in B){
			r<-regexec(pattern="([0-9]*)\\s*[*]?\\s*(\\w+)",text=b)
			m <- unlist(regmatches(b,r))
			print(m)
			RefName <- m[3]
			spr  <-  Reaction_createProduct(reaction)
			SimpleSpeciesReference_setSpecies(spr, RefName)
			if (nchar(m[2])>0){
			    SpeciesReference_setStoichiometry(spr,as.numeric(m[2]))
			}
		}
	}
}

.sbtab.expression.to.sbml <- function(sbml,Expression,CompartmentName,Compound){
	all.uid <- .unit.id.from.string(Expression$Unit)
	num.species <- nrow(Expression)
	ExpressionName <- row.names(Expression)
	for (i in 1:num.species){
		print(ExpressionName[i])
		this.unit.id <- all.uid[i]
		message(sprintf("unit of expression %i (%s): «%s»",i,ExpressionName[i],this.unit.id))
		F <- Expression$Formula[i]
		if (grepl("^\\s*[+-]?[0-9]+\\s*$",F)){
			message("this expression is a constant, mapping it to a parameter")
			par <- Model_createParameter(sbml)
			Parameter_setId(par,ExpressionName[i]);
			Parameter_setName(par,ExpressionName[i]);
			Parameter_setValue(par, as.numeric(F));
			Parameter_setUnits(par, this.unit.id);
		} else {
			pat <- paste0("^",ExpressionName[i],"$")
			if (any(grepl(pat,Compound$Assignment))){
			    j <- grep(pat,Compound$Assignment)
			    stopifnot(length(j)==1)
			    VarName <- row.names(Compound)[j]
			    if (is.na(VarName)) {
			        print(ExpressionName[i])
			        print(pat)
			        print(j)
			        print(row.names(Compound)[j])
			        print(row.names(Compound))
			        stop("Assignment not understood correctly.")
			    }
			} else {
			    sp <- Model_createSpecies(sbml);
			    Species_setId(sp, ExpressionName[i]);
			    Species_setName(sp, ExpressionName[i]);
			    Species_setCompartment(sp, CompartmentName);
			    Species_setBoundaryCondition(sp, "true")
			    VarName <- ExpressionName[i]
			}
			astMath <- parseL3Formula(F);
			##print(F)

			rule <- Model_createAssignmentRule(sbml)
			Rule_setVariable(rule,VarName)
			Rule_setMath(rule,astMath)
		}
	}
}

.sbtab.output.to.sbml <- function(sbml,Output,CompName,OutputUnitID){
	all.uid <- .unit.id.from.string(Output$Unit)
	num.species <- nrow(Output)
	Name <- row.names(Output)
	for (i in 1:num.species){
		message(sprintf("output %i:",i))
		print(Name[i])
		sp <- Model_createSpecies(sbml)
		this.unit.id <- all.uid[i]
		message(sprintf("unit of output %i (%s): «%s»",i,Name[i],this.unit.id))
		Species_setId(sp, Name[i])

		Species_setCompartment(sp, CompName)
		Species_setUnits(sp,OutputUnitID)
		Species_setInitialConcentration(sp, 0.0)
		Species_setBoundaryCondition(sp,"true")
		Species_setName(sp, Name[i])
				## SBase_appendAnnotation(sp,"<annotation>Output</annotation>")

		F <- Output$Formula[i]
		astMath <- parseL3Formula(F)
				##print(astMath)
		message(sprintf("output %i has formula:",i))
		print(F)
		message(formulaToString(astMath))
		rule <- Model_createAssignmentRule(sbml)
		Rule_setVariable(rule,Name[i])
		Rule_setMath(rule,astMath)
	}
}

#' Create an SBML model from Tables
#'
#' This function creates a character vector that contains the
#' biological model in SBML form. The return value can be written into
#' a file to create an SBML document.
#'
#' The input is an SBtab model as data.frames
#'
#' @param H model header, a Model identifying string
#' @param Defaults A data.frame with default units
#' @param Constant Table of Constants for this model [optional]
#' @param Parameter Table of Parameters
#' @param Input Table of input parameters, treated the same as normal Parameters
#' @param Expression will be interpreted and transformed into somthing that SBML understands (optional)
#' @param Reaction Table of reactions
#' @param Compound Table of species (reacting compounds)
#' @param Output Table of output functions (a model of measurement/observables)
#' @param Comp Table of Compartments
#' @import libSBML
#' @return SBML model as a character vector
.make.sbml <- function(H,Defaults,Constant=NULL,Parameter,Input=NULL,Expression=NULL,Reaction,Compound,Output,Comp){
	Doc <- SBMLDocument(2,4); # (LEVEL,VERSION)
	sbml  <-  SBMLDocument_createModel(Doc);
	Model_setId(sbml,H);
	## default units
	## substance
	for (i in 1:nrow(Defaults)){
		message(sprintf("default unit «%s» (%s)",Defaults$ID[i],Defaults$Unit[i]))
		unit <- .interpret.unit.from.string(Defaults$Unit[i])
		.create.sbml.unit(sbml,unit,Defaults$ID[i])
	}

	u.units <- .unique.units(sbml,c(Expression$Unit,Compound$Unit,Parameter$Unit,Constant$Unit,Input$Unit,Output$Unit,Comp$Unit))
	CompName <- row.names(Comp)
	print(CompName)
	for (i in 1:nrow(Comp)){
		this.unit.id <- .unit.id.from.string(Comp$Unit[i])
		comp <- Model_createCompartment(sbml)
		Compartment_setId(comp,Comp$ID[i])
		Compartment_setName(comp,CompName[i])
		Compartment_setUnits(comp,this.unit.id)
		Compartment_setVolume(comp,as.numeric(Comp$Size[i]))
	}
	if(!is.null(Constant)){
		.sbtab.constant.to.sbml(sbml,Constant)
	}
	.sbtab.compound.to.sbml(sbml,Compound,Comp$ID[1],"substance")
	.sbtab.parameter.to.sbml(sbml,Parameter)
	if (!is.null(Input)){
		.sbtab.parameter.to.sbml(sbml,Input)
	}
	.sbtab.reaction.to.sbml(sbml,Reaction,Comp)
	if (!is.null(Expression)){
		.sbtab.expression.to.sbml(sbml,Expression,Comp$ID[1],Compound)
	}
	if (!is.null(Output)){
		.sbtab.output.to.sbml(sbml,Output,Comp$ID[1],"substance")
	}

	return(Doc)
}
