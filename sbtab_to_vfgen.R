## this function expects a model to be loaded from an ODS file using
## read.ods or read_ods from the readODS package
##
##
GetTableName <- function(sbtab){
    ## the table title has to be in either one of the columns of row 1
    N <- names(sbtab);
    lN <- length(N);
    pattern <- "TableName='([^']+)'";
    m <- regexec(pattern, sbtab[1,], useBytes=TRUE);
    match <- regmatches(sbtab[1,],m);
    table.name <- "";
    for (i in c(1:lN)){
        if (length(match[[i]])>1){
            table.name <- match[[i]][2]; # so the first experssion in parentheses
        } 
    }
    return(table.name);
}

GetDocumentName <- function(sbtab){
    N <- names(sbtab);
    lN <- length(N);
    pattern <- "Document='([^']+)'";
    m <- regexec(pattern, sbtab[1,], useBytes=TRUE);
    match <- regmatches(sbtab[1,],m);
    document.name <- "Model";
    for (i in c(1:lN)){
        if (length(match[[i]])>1){
            document.name <- match[[i]][2]; # so the first experssion in parentheses
        } 
    }
    return(document.name);
}

make.cnames <- function(Labels){
    Unique.Names <- gsub("[.]","_",make.names(trimws(Labels), unique = TRUE, allow_ = TRUE),useBytes=TRUE);
    return(Unique.Names);
}

GetConservationLaws <- function(N){
    M <- pracma::null(t(N));
    if (dim(M)[1] > 1){
        M <- t(pracma::rref(t(M)));
    } else {
        M <- M/max(M);
    }
    nr=M;
    count=0;
    f <- c(2,3,5,7);
    while (norm(nr-round(nr),type="F") > 1e-6 && count<length(f)){
        count <- count+1;
        message(sprintf("nullspace is not represented by integers. \nTo make the mass conservation more readable, we multiply them by %i and round.",f[count]));
        nr <- nr*f[count];
    }
    Laws=round(nr);
    return(Laws);
}

AppendAmounts <- function(S,Quantity,QuantityName,Separator){
    n <- length(QuantityName);
    aq <- abs(as.numeric(Quantity));
    s <- paste(aq,QuantityName,sep="*",collapse=Separator);
    S <- paste0(S,s);
    return(S);
}

GetLawText <- function(Laws,CompoundName,InitialValue){
    nLaws <- dim(Laws)[2];
    nC <- length(CompoundName);
    I <- 1:nC;
    ConLaw <- list();
    ConLaw[["Constant"]] <- vector(length=nLaws);
    ConLaw[["Eliminates"]] <- vector(length=nLaws);
    ConLaw[["Formula"]] <- vector(length=nLaws,mode="character");
    ConLaw[["ConstantName"]] <- vector(length=nLaws,mode="character");
    ConLaw[["Text"]] <- vector(length=nLaws,mode="character");
    a <- vector(length=nC,mode="logical");
    a[] <- TRUE; # allowed replacement index...everything that is TRUE may be replaced by a conservation law expression
    for (j in 1:nLaws){
        l <- as.integer(round(Laws[,nLaws-j+1]));
        p <- l>0;
        n <- l<0;
        LawTextP <- AppendAmounts("",l[p],CompoundName[p],"+");
        LawTextN <- AppendAmounts("",l[n],CompoundName[n],"-");
        ##message(sprintf("length of LawTextN: %i («%s»)\n",nchar(LawTextN),LawTextN));
        if (nzchar(LawTextN)){
            LawText <- paste(LawTextP,LawTextN,sep="-");
        }else{
            LawText <- LawTextP;
        }        
        ## MaximumInitialValue <- max(InitialValue[(p|n)&a]); # the thing to be replace has to appear in the law, but not replaced by previous laws.
        
        u <- which.max(InitialValue[(p|n)&a])
        k <- I[(p|n)&a][u];
        ##print(k)
        ##message(sprintf("Replacing compound %i («%s») by Conservation Law Expression.\n",k,CompoundName[k]));
        a[k] <- FALSE; # this compound may not be replaced in the future;
        Const <- as.numeric(l %*% InitialValue);
        ConLaw$ConstantName[j] <- sprintf("%s_ConservedConst",CompoundName[k]);
        ConLaw$Text[j] <- paste(ConLaw$ConstantName[j],sub("1[*]","",LawText),sep=" = ");
        m <- (1:nC != k);

        ConLaw$Constant[j] <- Const;
        ConLaw$Eliminates[j] <- k;

        FormulaP <- AppendAmounts("", l[p&m], CompoundName[p & m],"+");
        FormulaN <- AppendAmounts("", l[n&m], CompoundName[n & m],"-");
        if (nzchar(FormulaN)){
            Formula <- paste(FormulaP,FormulaN,sep="-")
        }else{
            Formula <- FormulaP;
        }
        ConLaw$Formula[j] <- sub("1[*]","",Formula);
        message(sprintf("This will comment out compund %i («%s», initial value: %g), Conserved Constant = %f\n",k,CompoundName[k],InitialValue[k],Const));
    }
    return(ConLaw)
}


sbtab_to_vfgen <- function(M){
    lM <- length(M);
    SBtab <- list();
    document.name <- GetDocumentName(M[[1]]);
    message(sprintf("Document Name: %s.",document.name));
    table.name <- list();
    for (i in c(1:lM)){
        table.name[[i]] <- GetTableName(M[[i]]);        
        ## table.title <- GetTableTitle(M[[i]]);
        SBtab[[i]] <- M[[i]][-c(1,2),];
        names(SBtab[[i]]) <- M[[i]][2,];        
    }
    cat(sprintf("SBtab has %i tables.\n",length(SBtab)));
    names(SBtab) <- table.name;
    message(table.name);
    message("The names of the SBtab list:");
    message(names(SBtab),sep=", ");
    message("The names of SBtab[[1]]:");
    message(colnames(SBtab[[1]]),sep=", ");
    message("assigning specific table columns to variables:");

    ReactionID <- SBtab[["Reaction"]][["!ID"]];
    ReactionFormula <- SBtab[["Reaction"]][["!ReactionFormula"]];
    ReactionName <- make.cnames(SBtab[["Reaction"]][["!Name"]])
    nR <- length(ReactionID);

    ConstantID <- SBtab[["Constant"]][["!ID"]];    
    ConstantName <- make.cnames(SBtab[["Constant"]][["!Name"]])
    ConstantValue <- SBtab[["Constant"]][["!Value"]];
    nConst <- length(ConstantID);

    ExpressionID <- SBtab[["Expression"]][["!ID"]];    
    ExpressionName <- make.cnames(SBtab[["Expression"]][["!Name"]])
    ExpressionFormula <- SBtab[["Expression"]][["!Formula"]];
    nExpression <- length(ExpressionID);

    
    CompoundID <- SBtab[["Compound"]][["!ID"]];
    CompoundName <- make.cnames(SBtab[["Compound"]][["!Name"]])
    InitialValue <- as.numeric(SBtab[["Compound"]][["!InitialValue"]]);
    names(InitialValue) <- CompoundName;
    nC <- length(CompoundID);
    message("Initial Values:");
    print(InitialValue);
    
    ParID <- SBtab[["Parameter"]][["!ID"]];
    ParName <- make.cnames(SBtab[["Parameter"]][["!Name"]])
    if (length(grep("!DefaultValue",colnames(SBtab[["Parameter"]])))>0){
        ParValue <- SBtab[["Parameter"]][["!DefaultValue"]];
    }else if (length(grep("!Value",colnames(SBtab[["Parameter"]])))>0){
        ParValue <- SBtab[["Parameter"]][["!Value"]];
    }else if (length(grep("!Mean",colnames(SBtab[["Parameter"]])))>0){
        ParValue <- SBtab[["Parameter"]][["!Mean"]];
    }else if (length(grep("!Median",colnames(SBtab[["Parameter"]])))>0){
        ParValue <- SBtab[["Parameter"]][["!Median"]];
    }
    nPar <- length(ParID);
    
    ## diregard parameters that have been previously auto generated, by Conservation of Mass analysis:
    

    FluxID <- SBtab[["Reaction"]][["!ID"]];
    Flux <- SBtab[["Reaction"]][["!KineticLaw"]];
    FluxName <- make.cnames(SBtab[["Reaction"]][["!Name"]]);
    nFlux <- length(FluxID);
    
    OutputID <- SBtab[["Output"]][["!ID"]];
    OutputName <- make.cnames(SBtab[["Output"]][["!Name"]]);
    OutputFormula <-  SBtab[["Output"]][["!Formula"]];
    nO <- length(OutputID);

    InputID <- SBtab[["Input"]][["!ID"]];
    InputName <- make.cnames(SBtab[["Input"]][["!Name"]]);
    InputDefaultValue <-  SBtab[["Input"]][["!DefaultValue"]];
    Disregard <- as.logical(SBtab[["Input"]][["!ConservationLaw"]]);
    if (any(is.na(Disregard))){
        Disregard <- as.logical(as.numeric(SBtab[["Input"]][["!ConservationLaw"]]));
    }    
    message("Some input parameters may be earlier detected Conservation Law constants: ");
    print(Disregard)
    if (length(Disregard)==length(InputID)){
        InputID <- InputID[!Disregard];
        InputName <- InputName[!Disregard];
        InputDefaultValue <- InputDefaultValue[!Disregard];
    }
    nInput <- length(InputID);

    
    lhs_rhs <- strsplit(ReactionFormula,"<=>");
    ODE<-vector(mode="character",length=nC);
    ## stoichiometry matrix:
    N <- matrix(0,nrow=nC,ncol=nFlux);
    
    for (i in c(1:nR)){
        line=lhs_rhs[[i]];
        message(sprintf("Reaction %i:",i));
        cat(sprintf("line (a->b): «%s» → «%s»\n",line[1],line[2]));
        a <- unlist(strsplit(line[1],"[+]"));
        b <- unlist(strsplit(line[2],"[+]"));
        message(" where a: ");
        print(a)
        message("   and b: ");
        print(b)
        message("Products:")        
        for (xcompound in b){
            ## message(xcompound)
            ## find possible factors within string
            xb <- unlist(strsplit(trimws(xcompound),"[* ]"));
            if (length(xb)>1){
                n <- round(as.numeric(xb[1]));
                compound <- make.cnames(xb[2]);
            } else {
                compound <- make.cnames(xb[1]);
                n <- 1;
            }
            cat(sprintf("%i × %s",n,compound))
            if (compound %in% CompoundName){
                j <- match(compound,CompoundName)
                message(sprintf("\t\t\t(%s is compound %i)",compound,j));
                ODE[j] <- paste(ODE[j],FluxName[i],sep="+");
                N[j,i] <- N[j,i] + n;
            } else if (compound %in% ExpressionName){
                message(sprintf("\t\t\t%s is a fixed expression, it has no influx. ODE will be unaffected, but the expression may be used in ReactionFlux calculations\n",compound));
            } else if (compound %in% "null") {
                message(sprintf("\t\t\t%s (Ø) is a placeholder to formulate degradation in reaction formulae.\n",compound));
            } else {
                stop(sprintf("\t\t\t%s is neither in the list of registered compounds nor is it an expression\n",compound));                
            }
        }
        message("Reactants:")
        for (xcompound in a){
            message(compound)
            xa <- unlist(strsplit(trimws(xcompound),"[* ]"));
            if (length(xa)>1){
                n <- round(as.numeric(xa[1]));
                compound <- make.cnames(xa[2]);
            } else {
                compound <- make.cnames(xa[1]);
                n <- 1;
            }
            cat(sprintf("%i × %s",n,compound))

            if (compound %in% CompoundName){
                j <- match(compound,CompoundName)
                message(sprintf("\t\t\t(%s is compound %i)",compound,j));
                ODE[j] <- paste(ODE[j],FluxName[i],sep="-");
                N[j,i] <- N[j,i] - n;
            } else if (compound %in% "null") {
                message(sprintf("\t\t\t%s (Ø) is a placeholder to formulate degradation in reaction formulae.\n",compound));
            } else if (compound %in% ExpressionName){
                message(sprintf("\t\t\t%s is a fixed expression, it has no influx. ODE will be unaffected, but the expression may be used in ReactionFlux calculations\n",compound));
            } else {
                stop(sprintf("\t\t\t%s is neither in the list of registered compounds nor is it an expression\n",compound));                
            }
        }
    }

    Laws <- GetConservationLaws(N);
    nLaws <- dim(Laws)[2];
    
    message(sprintf("Number of compunds:\t%i\nNumber of Reactions:\t%i",nC,nFlux));
    message(sprintf("Conservation Law dimensions:\t%i × %i\n",dim(Laws)[1],dim(Laws)[2]));
    message(sprintf("To check that the conservation laws apply: norm(t(StoichiometryMatrix) * ConservationLaw == %6.5f)",norm(t(N) %*% Laws),type="F"));

    ## currently the laws are not used, but each line of Laws will replace one state variable.
    print(t(Laws))
    
    ConLaw <- GetLawText(Laws,CompoundName,InitialValue);
    message(sprintf("!ID;!Name;!DefaultValue;!Unit;!ConservationLaw;!Formula"))
    for (j in 1:nLaws){
        message(sprintf("CLU%i;%s;%g;nM;TRUE;%s",j,ConLaw$ConstantName[j],ConLaw$Constant[j],ConLaw$Text[j]));
    }
    message(sprintf("If you'd like to monitor omitted compounds, add this to the Output table:\n"));
    message(sprintf("!ID;!Name;!Comment;!ErrorName;!ErrorType;!Unit;!ProbDist;!Formula"))
    for (j in 1:nLaws){
        k <- ConLaw$Eliminates[j];
        message(sprintf("YCL%i;%s_mon;monitors implicit state;SD_YCL%i;not applicable;nM;none;%s",j,CompoundName[k],j,CompoundName[k]));
    }
    
    H <- document.name;
    H <- sub("-",'_',H);
    vfgen.header <- "<?xml version=\"1.0\" ?>";
    vfgen.model <- sprintf("<VectorField Name=\"%s\" Description=\"model created by an R script «sbtab_to_vfgen.R»\">",H);
    vfgen.const <- vector(length=nConst,mode="character");
    for (i in c(1:nConst)){
        vfgen.const[i] <- sprintf(" <Constant Name=\"%s\" Description=\"constant %s\" Value=\"%s\"/>",ConstantName[i],ConstantID[i],ConstantValue[i]);
    }
    vfgen.par <- vector(length=nPar,mode="character");
    for (i in c(1:nPar)){
        vfgen.par[i] <- sprintf(" <Parameter Name=\"%s\" Description=\"independent parameter %s\" DefaultValue=\"%s\"/>",ParName[i],ParID[i],ParValue[i]);
    }
    vfgen.input <- vector(length=nInput,mode="character");
    for (i in c(1:nInput)){
        vfgen.input[i] <- sprintf(" <Parameter Name=\"%s\" Description=\"input parameter %s\" DefaultValue=\"%s\"/>",InputName[i],InputID[i],InputDefaultValue[i]);
    }

    vfgen.ConservationInput <- vector(length=nLaws,mode="character");
    for (i in c(1:nLaws)){
        k <- ConLaw$Eliminates[i];        
        vfgen.ConservationInput[i] <- sprintf(" <Parameter Name=\"%s\" Description=\"conserved quantity; eliminates %s as a state variable\" DefaultValue=\"%f\"/>",ConLaw$ConstantName[i],CompoundName[k],ConLaw$Constant[i]);
    }
    vfgen.ConservationLaw <- vector(length=nLaws,mode="character");
    for (i in c(1:nLaws)){
        k <- ConLaw$Eliminates[i];
        ##print(k)
        ##message(sprintf("Replacing compound %i («%s») by Conservation Law Expression.\n",k,CompoundName[k]));        
        F <- sprintf("%s - (%s)",ConLaw$ConstantName[i],ConLaw$Formula[i]);
        vfgen.ConservationLaw[i] <- sprintf(" <Expression Name=\"%s\" Description=\"derived from conservation law %i\" Formula=\"%s\"/>",CompoundName[k],i,F);
    }
    vfgen.expression=vector(length=nExpression,mode="character");
    for (i in c(1:nExpression)){
        vfgen.expression[i] <- sprintf(" <Expression Name=\"%s\" Description=\"defined expression %s\" Formula=\"%s\"/>",ExpressionName[i],ExpressionID[i],ExpressionFormula[i]);
    }    
    vfgen.flux=vector(length=nFlux,mode="character");
    for (i in c(1:nFlux)){
        vfgen.flux[i] <- sprintf(" <Expression Name=\"%s\" Description=\"flux %s\" Formula=\"%s\"/>",FluxName[i],FluxID[i],Flux[i]);
    }
    vfgen.ode=vector(length=nC,mode="character");
    for (i in c(1:nC)){
        if (i %in% ConLaw$Eliminates){
            message(sprintf("StateVariable %s will be commented out as it was laready defined as a Mass Conservation Law Expression.",CompoundName[i]));
            vfgen.ode[i] <- sprintf("<!-- <StateVariable Name=\"%s\" Description=\"compound %s\" DefaultInitialCondition=\"%s\" Formula=\"%s\"/> -->",CompoundName[i], CompoundID[i], InitialValue[i],ODE[i]);
        }else{
            vfgen.ode[i] <- sprintf(" <StateVariable Name=\"%s\" Description=\"compound %s\" DefaultInitialCondition=\"%s\" Formula=\"%s\"/>",CompoundName[i], CompoundID[i], InitialValue[i],ODE[i]);
        }
    }
    vfgen.function=vector(length=nO,mode="character");
    for (i in c(1:nO)){
        vfgen.function[i]=sprintf(" <Function Name=\"%s\" Description=\"output %s\" Formula=\"%s\"/>",OutputName[i],OutputID[i],OutputFormula[i]);
    }
    vfgen.endmodel <- "</VectorField>";
    fname<-sprintf("%s_vf.xml",H);
    cat(vfgen.header,vfgen.model,vfgen.const,vfgen.par,vfgen.input,vfgen.ConservationInput,vfgen.ConservationLaw,vfgen.expression,vfgen.flux,vfgen.ode,vfgen.function,vfgen.endmodel,sep="\n",file=fname);

    message(sprintf("The content was written to: %s\n",fname));
    return(c(vfgen.header,vfgen.model,vfgen.const,vfgen.par,vfgen.input,vfgen.ConservationInput,vfgen.ConservationLaw,vfgen.expression,vfgen.flux,vfgen.ode,vfgen.function,vfgen.endmodel));
}
