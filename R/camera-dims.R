#' @title Modified version of CAMERA so that it can be run on DI-MS dataset
#'
#' @description
#' Perform the CAMERA annotation functions on a an DI-MS dataset consisting of a column of mz values and a
#' corresponding column of intensity values.
#'
#' This was a quick hack as we wanted to get CAMERA working on DI-MS without having to make an official change to the software.
#'
#' Full documentation for CAMERA can be found here
#' https://bioconductor.org/packages/release/bioc/html/CAMERA.html
#'
#' @param data data.frame = two column data frame consisting of mz values and intensity
#' @param params_iso list = a list of all the parameters used for the isotope calculations
#' @param params_iso list = a list of all the parameters used for the adduct calculations
#' @param rule_type character = [primary, extended, user]
#' @param rule_pth character = path of the rules used for adducts. only required for rule type 'user'
#' @param rule_sep character = seperator for rulepath
#' @param rule_export boolean = If TRUE will export the ruleset used
#' @export
#'
cameraDIMS <- function(data, params_iso, params_adduct, rule_type='extended', rule_pth=NULL, rule_sep='\t',
                       rule_export=FALSE){
  if(nrow(data.frame(data))<=1){
    data <- data.frame(data, row.names=NULL)
    data$isotopes = ""
    data$adduct = ""
    data$pcgroup = 1
    return(list("peaklist"=data, "annoID"=NA))
  }

  if(params_adduct$polarity=="pos"){
    params_adduct$polarity = "+"
  }else if (params_adduct$polarity=="neg") {
    params_adduct$polarity = "-"
  }

  # remove any previous annotation
  data <- data[ , -which(names(data) %in% c("isotopes","adduct"))]


  if(rule_type=='extended'){

    if(params_adduct$polarity=="+"){
      ruleF <- read.table(system.file("rules", "extended_adducts_pos.txt", package = "cameraDIMS"), sep=rule_sep,
                          header=TRUE)
    }else if (params_adduct$polarity=="-") {
      ruleF <- read.table(system.file("rules", "extended_adducts_neg.txt", package = "cameraDIMS"), sep=rule_sep,
                          header=TRUE)
    }
  }else if (rule_type=='primary'){
    print('CHECKKKK')
    print(params_adduct$polarity)
    if(params_adduct$polarity=="+"){
      ruleF <- read.table(system.file("rules", "primary_adducts_pos.txt", package = "cameraDIMS"), sep=rule_sep,
                          header=TRUE)
    }else if (params_adduct$polarity=="-") {
      ruleF <- read.table(system.file("rules", "primary_adducts_neg.txt", package = "cameraDIMS"), sep=rule_sep,
                          header=TRUE)
    }
  }else if (rule_type=='user'){
    if(!is.null(rule_pth)){
      ruleF <- read.table(rule_pth, header = TRUE, sep=rule_sep)
    }else{
      print('rule_type == user, then a valid rule_pth is required')
      return(0)
    }
  }else{
    print('rule_type needs to be either [primary, extended, user]')
    return(0)
  }


  print('number of adduct rules:')
  print(nrow(ruleF))
  print('rules head:')
  print(head(ruleF))
  data <- data[order(data$mz),]
  print('head data:')
  print(head(data))

  # Do the isotope annotation
  isotopes <- dims_isotopes(data, params_iso)

  # do the adduct annotation
  annotated_out <- dims_adducts(data, ruleF, isotopes, params_adduct)

  if (rule_export){
    annotated_out[[3]] = ruleF
  }

  return(annotated_out)
}


dims_isotopes <- function(data, params){
  print("Performing isotope annotation")
  if(is.null(params$IM)){
    params$IM <- calcIsotopeMatrix()
  }else{
    isotopeMatrix <- params$IM
  }

  #number of peaks in pseudospectra
  ncl <- nrow(data)
  imz <- data$mz
  irt = rep(1, nrow(data))
  mint <- data$i
  isotope <- vector("list", length(imz));
  isomatrix <- matrix(ncol=5, nrow=0);
  colnames(isomatrix) <- c("mpeak", "isopeak", "iso", "charge", "intrinsic")
  mz <- imz
  int <- matrix(mint)
  ipeak <- seq(1, length(mz))

  #########################################
  # Find isotope for fake pseudo spectrum
  #########################################
  isomatrix <- findIsotopesPspec(isomatrix, mz, ipeak, int, params)

  #########################################
  # Run remaining findIsotope camera function
  #########################################
  #clean isotopes
  if(is.null(nrow(isomatrix))) {
    isomatrix = matrix(isomatrix, byrow=F, ncol=length(isomatrix))
  }
  #check if every isotope has only one annotation
  if(length(idx.duplicated <- which(duplicated(isomatrix[, 2]))) > 0){
    peak.idx <- unique(isomatrix[idx.duplicated, 2]);
    for( i in 1:length(peak.idx)){
      #peak.idx has two or more annotated charge
      #select the charge with the higher cardinality
      peak <- peak.idx[i];
      peak.mono.idx <- which(isomatrix[,2] == peak)
      if(length(peak.mono.idx) < 2){
        #peak has already been deleted
        next;
      }
      peak.mono <- isomatrix[peak.mono.idx,1]
      #which charges we have
      charges.list <- isomatrix[peak.mono.idx, 4];
      tmp <- cbind(peak.mono,charges.list);
      charges.length <- apply(tmp,1, function(x,isomatrix) {
        length(which(isomatrix[, 1] == x[1] & isomatrix[,4] == x[2])) },
        isomatrix);
      idx <- which(charges.length == max(charges.length));
      if(length(idx) == 1){
        #max is unique
        isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }else{
        #select this one, which lower charge
        idx <- which.min(charges.list[idx]);
        isomatrix <- isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }
    }
  }


  #check if every isotope in one isotope grp, have the same charge
  if(length(idx.duplicated <- which(duplicated(paste(isomatrix[, 1], isomatrix[, 3])))) > 0){
    #at least one pair of peakindex and number of isotopic peak is identical
    peak.idx <- unique(isomatrix[idx.duplicated,1]);
    for( i in 1:length(peak.idx)){
      #peak.idx has two or more annotated charge
      #select the charge with the higher cardinality
      peak <- peak.idx[i];
      #which charges we have
      charges.list <- unique(isomatrix[which(isomatrix[, 1] == peak), 4]);
      #how many isotopes have been found, which this charges
      charges.length <- sapply(charges.list, function(x,isomatrix,peak) { length(which(isomatrix[, 1] == peak & isomatrix[, 4] == x)) },isomatrix,peak);
      #select the charge which the highest cardinality
      idx <- which(charges.length == max(charges.length));
      if(length(idx) == 1){
        #max is unique
        isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }else{
        #select this one, which lower charge
        idx <- which.min(charges.list[idx]);
        isomatrix <- isomatrix[-which(isomatrix[, 1] == peak & isomatrix[, 4] %in% charges.list[-idx]),, drop=FALSE]
      }
    }
  }
  #Combine isotope cluster, if they overlap
  index2remove <- c();
  if(length(idx.duplicated <- which(isomatrix[, 1] %in% isomatrix[, 2]))>0){
    for(i in 1:length(idx.duplicated)){
      index <- which(isomatrix[, 2] == isomatrix[idx.duplicated[i], 1])
      index2 <- sapply(index, function(x, isomatrix) which(isomatrix[, 1] == isomatrix[x, 1] & isomatrix[,3] == 1),isomatrix)
      if(length(index2) == 0){
        index2remove <- c(index2remove,idx.duplicated[i])
      }
      max.index <- which.max(isomatrix[index,4]);
      isomatrix[idx.duplicated[i], 1] <- isomatrix[index[max.index], 1];
      isomatrix[idx.duplicated[i], 3] <- isomatrix[index[max.index], 3]+1;
    }
  }
  if(length(index <- which(isomatrix[,"iso"] > params$maxiso)) > 0){
    index2remove <- c(index2remove, index)
  }
  if(length(index2remove) > 0){
    isomatrix <- isomatrix[-index2remove,, drop=FALSE];
  }


  isomatrix <- isomatrix[order(isomatrix[,1]),,drop=FALSE]
  #Create isotope matrix within object
  #object@isoID <- matrix(nrow=0, ncol=4);
  #colnames(object@isoID) <- c("mpeak", "isopeak", "iso", "charge");
  #Add isomatrix to object
  #object@isoID <- rbind(object@isoID, isomatrix[, 1:4]);
  # counter for isotope groups
  globalcnt <- 0;
  oldnum <- 0;
  if(nrow(isomatrix) > 0){
    for( i in 1:nrow(isomatrix)){
      if(!isomatrix[i, 1] == oldnum){
        globalcnt <- globalcnt+1;
        isotope[[isomatrix[i, 1]]] <- list(y=globalcnt, iso=0, charge=isomatrix[i, 4], val=isomatrix[i, 5]);
        oldnum <- isomatrix[i, 1];
      };
      isotope[[isomatrix[i,2]]] <- list(y=globalcnt,iso=isomatrix[i,3],charge=isomatrix[i,4],val=isomatrix[i,5]);
    }
  }
  return(isotope)

}


dims_adducts <- function(data, ruleF, isotopes, params_adduct){
  print("Performing adduct annotation")
  #default polarity set to positive
  polarity <- params_adduct$polarity

  # get parameters
  maxCharge <- params_adduct$maxCharge
  maxMol <- params_adduct$maxMol
  quasimolion <- params_adduct$quasimolion
  devppm <- params_adduct$devppm
  mzabs <- params_adduct$mzabs
  imz <- data$mz
  ipeak <- matrix(seq(1, nrow(data)))

  # Indicate with NA where isotopes have already been called
  for(x in seq(along = isotopes)){
    if(!is.null(isotopes[[x]])){
      if(isotopes[[x]]$iso != 0){
        imz[x] <- NA;
      }
    }
  }

  # Get the rules we want to use
  #rules.idx = row.names(subset(ruleF, nmol <= maxMol & charge <= maxCharge & quasi <= 1))
  rules.idx = row.names(ruleF[(ruleF$nmol <= maxMol) & (ruleF$charge <= maxCharge), ])
  ruleF = ruleF[rules.idx,]

  # get hypotheses of adducts
  print("--Getting annotation hypothesis")
  hypothese <- annotateGrp(ipeak, imz, ruleF, mzabs, devppm, isotopes, quasimolion, rules.idx=rules.idx)

  if(is.null(hypothese)){
    #nothing found
    if(is.null(unlist(isotopes))){
      isotopesv <- ""
    }else{
      isotopesv <- vector("character", nrow(data));
      for(i in 1:length(isotopes)){
        if(length(isotopes) > 0&& !(is.null(isotopes[[i]]))) {
          num.iso <- isotopes[[i]]$iso;
          #Which isotope peak is peak i?
          if(num.iso == 0){
            str.iso <- "[M]";
          } else {
            str.iso <- paste("[M+", num.iso, "]", sep="")
          }
          #Multiple charged?
          if(isotopes[[i]]$charge > 1){
            isotopesv[i] <- paste("[", isotopes[[i]]$y, "]", str.iso, isotopes[[i]]$charge, polarity, sep="");
          }else{
            isotopesv[i] <- paste("[", isotopes[[i]]$y, "]", str.iso, polarity, sep="");
          }
        } else {
          #No isotope informationen available
          isotopesv[i] <- "";
        }
      }
    }

    rownames(data)<-NULL;#Bugfix for: In data.row.names(row.names, rowsi, i) :  some row.names duplicated:
    final_table <- cbind(data, isotopes = isotopesv, adduct="", pcgroup=1)

    return(list("peaklist"=final_table, "annoID"=NA))


  }

  oidscore <- c();
  index <- c();
  annoID <- matrix(ncol=4, nrow=0)
  annoGrp <- matrix(ncol=4, nrow=0)

  # I think we can get 'parent' to work by simply changing the annoID group column name from parentID to parent
  #colnames(annoID) <-  c("id", "grpID", "ruleID", "parent")

  colnames(annoID) <-  c("id", "grpID", "ruleID", "parentID")
  colnames(annoGrp) <- c("id", "mass", "ips", "psgrp")
  charge <- 0;
  massgrp <- 0;
  old_massgrp <- 0
  #old_massgrp <- min(hypothese[,'massgrp']);
  parent <- FALSE
  i <- 1 #pspectra_list[j]

  for(hyp in 1:nrow(hypothese)){
    peakid <- as.numeric(ipeak[hypothese[hyp, "massID"]]);
    if(old_massgrp != hypothese[hyp, "massgrp"]) {
      massgrp <- massgrp + 1;
      old_massgrp <- hypothese[hyp, "massgrp"];
      annoGrp <- rbind(annoGrp, c(massgrp, hypothese[hyp, "mass"],
                                  sum(hypothese[ which(hypothese[, "massgrp"] == old_massgrp), "score"]), i) )
    }
    if(parent){
      annoID <- rbind(annoID, c(peakid, massgrp, hypothese[hyp, "ruleID"], ipeak[hypothese[hyp, "parent"]]))
    }else{
      annoID <- rbind(annoID, c(peakid, massgrp, hypothese[hyp, "ruleID"], NA))
    }
  }


  print("--Get derivative ions")
  derivativeIons <- getderivativeIons(annoID, annoGrp, ruleF, length(imz));

  # Replace the isotopes taken as part of the object e.g. object@isotopes is now isotopes_obj
  # Prevents too many uses of the name isotopes
  isotopes_obj <- isotopes

  #allocate variables for CAMERA output
  adduct <- vector("character", nrow(data));
  isotopes <- vector("character", nrow(data));
  pcgroup <- vector("character", nrow(data));


  print("--join isotope and adduct information")
  #First isotope informationen and adduct informationen
  for(i in seq(along = isotopes)){
    #check if adduct annotation is present for peak i
    if(length(derivativeIons) > 0 && !(is.null(derivativeIons[[i]]))) {
      #Check if we have more than one annotation for peak i
      if(length(derivativeIons[[i]]) > 1) {
        #combine ion species name and rounded mass hypophysis
        names <- paste(derivativeIons[[i]][[1]]$name, signif(derivativeIons[[i]][[1]]$mass, 12));
        for(ii in 2:length(derivativeIons[[i]])) {
          names <- paste(names, derivativeIons[[i]][[ii]]$name, signif(derivativeIons[[i]][[ii]]$mass, 12));
        }
        #save name in vector adduct
        adduct[i] <- names;
      } else {
        #Only one annotation
        adduct[i] <- paste(derivativeIons[[i]][[1]]$name, signif(derivativeIons[[i]][[1]]$mass, 12));
      }
    } else {
      #no annotation empty name
      adduct[i] <- "";
    }
    #Check if we have isotope informationen about peak i
    if(length(isotopes_obj) > 0&& !(is.null(isotopes_obj[[i]]))) {
      num.iso <- isotopes_obj[[i]]$iso;
      #Which isotope peak is peak i?
      if(num.iso == 0){
        str.iso <- "[M]";
      } else {
        str.iso <- paste("[M+", num.iso, "]", sep="")
      }
      #Multiple charged?
      if(isotopes_obj[[i]]$charge > 1){
        isotopes[i] <- paste("[", isotopes_obj[[i]]$y, "]", str.iso, isotopes_obj[[i]]$charge, polarity, sep="");
      }else{
        isotopes[i] <- paste("[", isotopes_obj[[i]]$y, "]", str.iso, polarity, sep="");
      }
    } else {
      #No isotope informationen available
      isotopes[i] <- "";
    }
  }
  #Have we more than one pseudospectrum?
  pcgroup <- 1

  rownames(data)<-NULL;#Bugfix for: In data.row.names(row.names, rowsi, i) :  some row.names duplicated:
  final_table <- data.frame(data, isotopes, adduct, pcgroup, stringsAsFactors=FALSE, row.names=NULL)
  return(list("peaklist"=final_table, "annoID"=annoID))
}
