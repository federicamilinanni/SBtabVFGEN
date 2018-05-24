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
    M <- MASS::Null(t(N));
    rM <- t(pracma::rref(t(M)));
    nr=rM;
    count=0;
    while (norm(nr-round(nr)) > 1e-6 && count<3){
        warning("nullspace is not represented by integers. \nTo make the mass conservation more readable, we multiply them by 10 and round.");
        nr <- nr*10;
        count <- count+1;
    }
    Laws=round(nr);
    return(Laws);
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

    CompoundID <- SBtab[["Compound"]][["!ID"]];
    CompoundName <- make.cnames(SBtab[["Compound"]][["!Name"]])
    InitialValue <- SBtab[["Compound"]][["!InitialValue"]];
    nC <- length(CompoundID);

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

    FluxID <- SBtab[["Reaction"]][["!ID"]];
    Flux <- SBtab[["Reaction"]][["!KineticLaw"]];
    FluxName <- make.cnames(SBtab[["Reaction"]][["!Name"]]);
    nFlux <- length(FluxID);
    
    OutputID <- SBtab[["Output"]][["!ID"]];
    OutputName <- make.cnames(SBtab[["Output"]][["!QuantityName"]]);
    OutputFormula <-  SBtab[["Output"]][["!Formula"]];
    nO <- length(OutputID);
    
    lhs_rhs <- strsplit(ReactionFormula,"<=>");
    ODE<-vector(mode="character",length=nC);
    ## stoichiometry matrix:
    N <- matrix(0,nrow=nC,ncol=nFlux);
    
    for (i in c(1:nR)){
        line=lhs_rhs[[i]];
        message(sprintf("Reaction %i:",i));
        cat(sprintf("line (a->b): «%s» -> «%s»\n",line[1],line[2]));
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
            xb <- unlist(strsplit(xcompound,"[*]"));
            if (length(xb)>1){
                n <- as.numeric(xb[1]);
                compound <- make.cnames(xb[2]);
            } else {
                compound <- make.cnames(xb[1]);
                n <- 1;
            }            
            message(sprintf("%i × %s",n,compound))
            if (compound %in% CompoundName){
                j <- match(compound,CompoundName)
                message(sprintf("%s is compound %i\n",compound,j));
                ODE[j] <- paste(ODE[j],FluxName[i],sep="+");
                N[j,i] <- N[j,i] + n;
            }
        }
        message("Reactants:")
        for (xcompound in a){
            message(compound)
            xa <- unlist(strsplit(xcompound,"[*]"));
            if (length(xa)>1){
                n <- as.numeric(xa[1]);
                compound <- make.cnames(xa[2]);
            } else {
                compound <- make.cnames(xa[1]);
                n <- 1;
            }
            message(sprintf("%i × %s",n,compound))

            if (compound %in% CompoundName){
                j <- match(compound,CompoundName)
                message(sprintf("%s is compound %i\n",compound,j));
                ODE[j] <- paste(ODE[j],FluxName[i],sep="-");
                N[j,i] <- N[j,i] - n;
            }
        }
    }

    Laws <- GetConservationLaws(N);
    print(Laws);
    ## currently the laws are not used, but each line of Laws will replace one state variable.
    
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
    vfgen.flux=vector(length=nFlux,mode="character");
    for (i in c(1:nFlux)){
        vfgen.flux[i] <- sprintf(" <Expression Name=\"%s\" Description=\"flux %s\" Formula=\"%s\"/>",FluxName[i],FluxID[i],Flux[i]);
    }
    vfgen.ode=vector(length=nC,mode="character");
    for (i in c(1:nC)){
        vfgen.ode[i] <- sprintf(" <StateVariable Name=\"%s\" Description=\"compound %s\" DefaultInitialCondition=\"%s\" Formula=\"%s\"/>",CompoundName[i], CompoundID[i], InitialValue[i],ODE[i]);
    }
    vfgen.function=vector(length=nO,mode="character");
    for (i in c(1:nO)){
        vfgen.function[i]=sprintf(" <Function Name=\"%s\" Description=\"output %s\" Formula=\"%s\"/>",OutputName[i],OutputID[i],OutputFormula[i]);
    }
    vfgen.endmodel <- "</VectorField>";
    fname<-sprintf("%s_vf.xml",H);
    cat(vfgen.header,vfgen.model,vfgen.const,vfgen.par,vfgen.flux,vfgen.ode,vfgen.function,vfgen.endmodel,sep="\n",file=fname);

    message(sprintf("The content was written to: %s\n",fname));
    return(c(vfgen.header,vfgen.model,vfgen.const,vfgen.par,vfgen.flux,vfgen.ode,vfgen.function,vfgen.endmodel));
}
