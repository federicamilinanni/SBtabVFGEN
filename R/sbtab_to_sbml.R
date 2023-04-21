.create.sbml.unit <- function(sbml,unit,id){
	n <- nrow(unit)
	unitdef <- Model_createUnitDefinition(sbml);
	UnitDefinition_setId(unitdef,id);
	for (i in 1:n){
		u <- UnitDefinition_createUnit(unitdef);
		Kind <- switch(unit$kind[i],litre="UNIT_KIND_LITRE",metre="UNIT_KIND_METRE",second="UNIT_KIND_SECOND",mole="UNIT_KIND_MOLE",kilogram="UNIT_KIND_KILOGRAM",gram="UNIT_KIND_GRAM","UNIT_KIND_DIMENSIONLESS")
		Unit_setKind(u, Kind);
		Unit_setExponent(u,unit$exponent[i]);
		Unit_setMultiplier(u,unit$multiplier[i]);
		Unit_setScale(u,unit$scale[i]);
	}
}

.unique.units <- function(sbml,Unit){
	U <- unique(Unit)
	uid <- unit.id(U,prnt=TRUE)
	unit <- list()
	for (i in 1:length(uid)){
		unit[[i]] <- unit.from.string(U[i])
		.create.sbml.unit(sbml,unit[[i]],uid[i])
	}
	names(unit) <- uid
	return(unit)
}

.sbtab.parameter.to.sbml <- function(sbml,Parameter){
	all.uid <- unit.id(Parameter$Unit)
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
	all.uid <- unit.id(Compound$Unit)
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
	all.uid <- unit.id(Constant$Unit)
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
	all.uid <- unit.id(Expression$Unit)
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
			rule <- Model_createAssignmentRule(sbml)
			Rule_setVariable(rule,VarName)
			Rule_setMath(rule,astMath)
		}
	}
}

.sbtab.output.to.sbml <- function(sbml,Output,CompName,OutputUnitID){
	all.uid <- unit.id(Output$Unit)
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
		F <- Output$Formula[i]
		astMath <- parseL3Formula(F)
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
		unit <- unit.from.string(Defaults$Unit[i])
		.create.sbml.unit(sbml,unit,Defaults$ID[i])
	}

	u.units <- .unique.units(sbml,c(Expression$Unit,Compound$Unit,Parameter$Unit,Constant$Unit,Input$Unit,Output$Unit,Comp$Unit))
	CompName <- row.names(Comp)
	print(CompName)
	for (i in 1:nrow(Comp)){
		this.unit.id <- unit.id(Comp$Unit[i])
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
