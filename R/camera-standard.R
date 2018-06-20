match_ranges <- function(IM, MI, max.index){
  match_m <- matrix(ncol=max.index, nrow = nrow(IM))
  for (i in 1:nrow(IM)){
    x <- IM[i,]
    ranges <- split(x, ceiling(seq_along(x)/2))
    df <- data.frame(t(data.frame(ranges)))
    colnames(df) <- c('start', 'end')
    df$id <- 1:nrow(df)


    out <- cut(MI[1:max.index], t(df[,c('start', 'end')]))
    levels(out)[c(FALSE,TRUE)]  <- NA
    match_m[i, ] <- df$id[out]


  }

  return(match_m)




}

calcIsotopeMatrix <- function(maxiso=4){
  if(!is.numeric(maxiso)){
    stop("Parameter maxiso is not numeric!\n")
  } else if(maxiso < 1 | maxiso > 8){
    stop(paste("Parameter maxiso must between 1 and 8. ",
               "Otherwise use your own IsotopeMatrix.\n"),sep="")
  }
  isotopeMatrix <- matrix(NA, 8, 4);
  colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")
  isotopeMatrix[1, ] <- c(1.000, 1.0040, 1.0, 150)
  isotopeMatrix[2, ] <- c(0.997, 1.0040, 0.01, 200)
  isotopeMatrix[3, ] <- c(1.000, 1.0040, 0.001, 200)
  isotopeMatrix[4, ] <- c(1.000, 1.0040, 0.0001, 200)
  isotopeMatrix[5, ] <- c(1.000, 1.0040, 0.00001, 200)
  isotopeMatrix[6, ] <- c(1.000, 1.0040, 0.000001, 200)
  isotopeMatrix[7, ] <- c(1.000, 1.0040, 0.0000001, 200)
  isotopeMatrix[8, ] <- c(1.000, 1.0040, 0.00000001, 200)
  return(isotopeMatrix[1:maxiso, , drop=FALSE])
}

findIsotopesPspec <- function(isomatrix, mz, ipeak, int, params){
  #isomatrix - isotope annotations (5 column matrix)
  #mz - m/z vector, contains all m/z values from specific pseudospectrum
  #int - int vector, see above
  #maxiso - how many isotopic peaks are allowed
  #maxcharge - maximum allowed charge
  #devppm - scaled ppm error
  #mzabs - absolut error in m/z
  #matrix with all important informationen
  spectra <- matrix(c(mz, ipeak), ncol=2)
  #int <- int[order(spectra[, 1]), , drop=FALSE]
  spectra <- spectra[order(spectra[, 1]), ];
  cnt <- nrow(spectra);
  #isomatrix <- matrix(NA, ncol=5, nrow=0)
  #colnames(isomatrix) <- c("mpeak", "isopeak", "iso", "charge", "intrinsic")
  #calculate error
  error.ppm <- params$devppm * mz;
  #error.abs <- ),1, function(x) x + params$mzabs*rbind(1,2,3)));
  #for every peak in pseudospectrum
  for ( j in 1:(length(mz) - 1)){
    #create distance matrix
    MI <- spectra[j:cnt, 1] - spectra[j, 1];
    #Sum up all possible/allowed isotope distances + error(ppm of peak mz and mzabs)
    max.index <- max(which(MI < (sum(params$IM[1:params$maxiso, "mzmax"]) + error.ppm[j] + params$mzabs )))
    #check if one peaks falls into isotope window
    if(max.index == 1){
      #no promising candidate found, move on
      next;
    }
    #IM - isotope matrix (column diffs(min,max) per charge, row num. isotope)
    IM <- t(sapply(1:params$maxcharge,function(x){
      mzmin <- (params$IM[, "mzmin"]) / x;
      mzmax <- (params$IM[, "mzmax"]) / x;
      error <- (error.ppm[j]+params$mzabs) / x
      res <- c(0,0);
      for(k in 1:length(mzmin)){
        res <- c(res, mzmin[k]+res[2*k-1], mzmax[k]+res[2*k])
      }
      res[seq(1,length(res),by=2)] <- res[seq(1,length(res),by=2)]-error
      res[seq(2,length(res),by=2)] <- res[seq(2,length(res),by=2)]+error
      return (res[-c(1:2)])
    } ))
    #Sort IM to fix bug, with high ppm and mzabs values
    #TODO: Find better solution and give feedback to user!
    IM <- t(apply(IM,1,sort))
    #find peaks, which m/z value is in isotope interval

    #hits <- t(apply(IM, 1, function(x){ findInterval(MI[1:max.index], x)}))
    hits <- match_ranges(IM, MI, max.index)

    rownames(hits) <- c(1:nrow(hits))
    colnames(hits) <- c(1:ncol(hits))
    hits[which(hits==0)] <-NA
    hits <- hits[, -1, drop=FALSE]
    #hits.iso <- hits%/%2 + 1;


    hits.iso <- hits


    #check occurence of first isotopic peak
    for(iso in 1:min(params$maxiso,ncol(hits.iso))){
      hit <- apply(hits.iso,1, function(x) any(naOmit(x)==iso))
      hit[which(is.na(hit))] <- TRUE
      if(all(hit)) break;
      hits.iso[!hit,] <- t(apply(hits.iso[!hit,,drop=FALSE],1, function(x) {
        if(!all(is.na(x))){
          ini <- which(x > iso)
          if(!is.infinite(ini) && length(ini) > 0){
            x[min(ini):ncol(hits.iso)] <- NA
          }
        }
        x
      }))
    }
    #set NA to 0
    hits[which(is.na(hits.iso))] <- 0
    #check if any isotope is found
    hit <- apply(hits, 1, function(x) sum(x)>0)
    #drop nonhits
    hits <- hits[hit, , drop=FALSE]
    #if no first isotopic peaks exists, next
    if(nrow(hits) == 0){
      next;
    }
    #getting max. isotope cluster length
    #TODO: unique or not????
    #isolength <- apply(hits, 1, function(x) length(which(unique(x) %% 2 !=0)))
    #isohits - for each charge, length of peak within intervals
    #isohits <- lapply(1:nrow(hits), function(x) which(hits[x, ] %% 2 !=0))
    isohits <- lapply(1:nrow(hits), function(x) which(hits[x, ] >0))
    isolength <- sapply(isohits, length)
    #Check if any result is found
    if(all(isolength==0)){
      next;
    }
    #itensity checks
    #candidate.matrix
    #first column - how often succeded the isotope intensity test
    #second column - how often could a isotope int test be performed
    candidate.matrix <- matrix(0, nrow=length(isohits), ncol=max(isolength)*2);
    for(iso in 1:length(isohits)){
      for(candidate in 1:length(isohits[[iso]])){
        for(sample.index in c(1:ncol(int))){
          #Test if C12 Peak is NA
          if(!is.na(int[j, sample.index])){
            #candidate.matrix[maxIso, 1] <- candidate.matrix[maxIso, 1] + 1
          }
          charge <- as.numeric(row.names(hits)[iso])
          int.c12 <- int[j, sample.index]
          isotopePeak <- hits[iso,isohits[[iso]][candidate]]%/%2 + 1;
          if(isotopePeak == 1){
            #first isotopic peak, check C13 rule
            int.c13 <- int[isohits[[iso]][candidate]+j, sample.index];
            int.available <- all(!is.na(c(int.c12, int.c13)))
            if (int.available){
              theo.mass <- spectra[j, 1] * charge; #theoretical mass
              numC <- abs(round(theo.mass / 12)); #max. number of C in molecule
              inten.max <- int.c12 * numC * 0.011; #highest possible intensity
              inten.min <- int.c12 * 1 * 0.011; #lowest possible intensity
              if((int.c13 < inten.max && int.c13 > inten.min) || !params$filter){
                candidate.matrix[iso,candidate * 2 - 1] <- candidate.matrix[iso,candidate * 2 - 1] + 1
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
              }else{
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
              }
            } else {
              #todo
            }
          } else {
            #x isotopic peak
            int.cx <- int[isohits[[iso]][candidate]+j, sample.index];
            int.available <- all(!is.na(c(int.c12, int.cx)))
            if (int.available) {
              intrange <- c((int.c12 * params$IM[isotopePeak,"intmin"]/100),
                            (int.c12 * params$IM[isotopePeak,"intmax"]/100))
              #filter Cx isotopic peaks muss be smaller than c12
              if(int.cx < intrange[2] && int.cx > intrange[1]){
                candidate.matrix[iso,candidate * 2 - 1] <- candidate.matrix[iso,candidate * 2 - 1] + 1
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
              }else{
                candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
              }
            } else {
              candidate.matrix[iso,candidate * 2 ] <- candidate.matrix[iso,candidate * 2] + 1
            }#end int.available
          }#end if first isotopic peak
        }#for loop samples
      }#for loop candidate
    }#for loop isohits
    #calculate ratios
    candidate.ratio <- candidate.matrix[, seq(from=1, to=ncol(candidate.matrix),
                                              by=2)] / candidate.matrix[, seq(from=2,
                                                                              to=ncol(candidate.matrix), by=2)];
    if(is.null(dim(candidate.ratio))){
      candidate.ratio <- matrix(candidate.ratio, nrow=nrow(candidate.matrix))
    }
    if(any(is.nan(candidate.ratio))){
      candidate.ratio[which(is.nan(candidate.ratio))] <- 0;
    }
    #decision between multiple charges or peaks
    for(charge in 1:nrow(candidate.matrix)){
      if(any(duplicated(hits[charge, isohits[[charge]]]))){
        #One isotope peaks has more than one candidate
        ##check if problem is still consistent
        for(iso in unique(hits[charge, isohits[[charge]]])){
          if(length(index <- which(hits[charge, isohits[[charge]]]==iso))== 1){
            #now duplicates next
            next;
          }else{
            #find best
            index2 <- which.max(candidate.ratio[charge, index]);
            save.ratio <- candidate.ratio[charge, index[index2]]
            candidate.ratio[charge,index] <- 0
            candidate.ratio[charge,index[index2]] <- save.ratio
            index <- index[-index2]
            isohits[[charge]] <- isohits[[charge]][-index]
          }
        }
      }#end if
      for(isotope in 1:ncol(candidate.ratio)){
        if(candidate.ratio[charge, isotope] >= params$minfrac){
          isomatrix <- rbind(isomatrix,
                             c(spectra[j, 2],
                               spectra[isohits[[charge]][isotope]+j, 2],
                               isotope, as.numeric(row.names(hits)[charge]), 0))
        } else{
          break;
        }
      }
    }#end for charge
  }#end for j
  return(isomatrix)
}
getderivativeIons <- function(annoID, annoGrp, rules, npeaks){
  #generate Vector length npeaks
  derivativeIons <- vector("list", npeaks);
  #intrinsic charge
  #TODO: Not working at the moment
  charge <- 0;
  #check if we have annotations
  if(nrow(annoID) < 1){
    return(derivativeIons);
  }
  for(i in 1:nrow(annoID)){
    peakid <- annoID[i, 1];
    grpid <- annoID[i, 2];
    ruleid <- annoID[i, 3];
    if(is.null(derivativeIons[[peakid]])){
      #Peak has no annotations so far
      if(charge == 0 | rules[ruleid, "charge"] == charge){
        derivativeIons[[peakid]][[1]] <- list( rule_id = ruleid,
                                               charge = rules[ruleid, "charge"],
                                               nmol = rules[ruleid, "nmol"],
                                               name = paste(rules[ruleid, "name"]),
                                               mass = annoGrp[grpid, 2])
      }
    } else {
      #Peak has already an annotation
      if(charge == 0 | rules[ruleid, "charge"] == charge){
        derivativeIons[[peakid]][[(length(
          derivativeIons[[peakid]])+1)]] <- list( rule_id = ruleid,
                                                  charge = rules[ruleid, "charge"],
                                                  nmol = rules[ruleid, "nmol"],
                                                  name=paste(rules[ruleid, "name"]),
                                                  mass=annoGrp[grpid, 2])
      }
    }
    charge <- 0;
  }
  return(derivativeIons);
}
getIsotopeCluster <- function(object, number=NULL, value="maxo",
                              sampleIndex=NULL){
  #check values
  if(is.null(object)) {
    stop("No xsa argument was given.\n");
  }else if(!class(object)=="xsAnnotate"){
    stop("Object parameter is no xsAnnotate object.\n");
  }
  value <- match.arg(value, c("maxo", "into", "intb"), several.ok=FALSE)
  if(!is.null(number) & !is.numeric(number)){
    stop("Number must be NULL or numeric");
  }
  if(!is.null(sampleIndex) & !all(is.numeric(sampleIndex))){
    stop("Parameter sampleIndex must be NULL or numeric");
  }
  if(is.null(sampleIndex)){
    nSamples <- 1;
  } else if( all(sampleIndex <= length(object@xcmsSet@filepaths) & sampleIndex > 0)){
    nSamples <- length(sampleIndex);
  } else {
    stop("All values in parameter sampleIndex must be lower equal
         the number of samples and greater than 0.\n")
  }
  if(length(sampnames(object@xcmsSet)) > 1){ ## more than one sample
    gvals <- groupval(object@xcmsSet, value=value);
    groupmat <- object@groupInfo;
    iso.matrix <- matrix(0, ncol=nSamples, nrow=length(object@isotopes));
    if(is.null(sampleIndex)){
      for(i in 1:length(object@pspectra)){
        iso.matrix[object@pspectra[[i]],1] <- gvals[object@pspectra[[i]],object@psSamples[i]];
      }
    } else {
      for(i in 1:length(object@pspectra)){
        iso.matrix[object@pspectra[[i]], ] <- gvals[object@pspectra[[i]], sampleIndex]
      }
    }
    peakmat <- cbind(groupmat[, "mz"], iso.matrix );
    rownames(peakmat) <- NULL;
    if(is.null(sampleIndex)){
      colnames(peakmat) <- c("mz",value);
    }else{
      colnames(peakmat) <- c("mz", sampnames(object@xcmsSet)[sampleIndex]);
    }
    if(any(is.na(peakmat))){
      cat("Warning: peak table contains NA values. To remove apply fillpeaks on xcmsSet.\n");
    }
  } else if(length(sampnames(object@xcmsSet)) == 1){ ## only one sample was
    peakmat <- object@groupInfo[, c("mz", value)];
  } else {
    stop("sampnames could not extracted from the xcmsSet.\n");
  }
  #collect isotopes
  index <- which(!sapply(object@isotopes, is.null));
  tmp.Matrix <- cbind(index, matrix(unlist(object@isotopes[index]), ncol=4, byrow=TRUE))
  colnames(tmp.Matrix) <- c("Index","IsoCluster","Type","Charge","Val")
  max.cluster <- max(tmp.Matrix[,"IsoCluster"])
  max.type <- max(tmp.Matrix[,"Type"])
  isotope.Matrix <- matrix(NA, nrow=max.cluster, ncol=(max.type+2));
  invisible(apply(tmp.Matrix,1, function(x) {
    isotope.Matrix[x["IsoCluster"],x["Type"]+2] <<- x["Index"];
    isotope.Matrix[x["IsoCluster"],1] <<- x["Charge"];
  }))
  invisible(apply(isotope.Matrix,1, function(x) {
    list(peaks=peakmat[na.omit(x[-1]),],charge=x[1])
  }))
  }


##Functions for findAdducts

annotateGrp <- function(ipeak, imz, rules, mzabs, devppm, isotopes, quasimolion, rules.idx) {
  #m/z vector for group i with peakindex ipeak
  mz     <- imz[ipeak];
  naIdx <- which(!is.na(mz))

  #Spectrum have only annotated isotope peaks, without monoisotopic peak
  #Give error or warning?
  if(length(na.omit(mz[naIdx])) < 1){
    return(NULL);
  }

  print("---Creating mass diff matrix")
  ML <- massDiffMatrix(mz[naIdx], rules[rules.idx,]);

  print("---createing hypothese")
  hypothese <- createHypothese(ML, rules[rules.idx, ], devppm, mzabs, naIdx);

  #create hypotheses
  if(is.null(nrow(hypothese)) || nrow(hypothese) < 2 ){
    return(NULL);
  }

  #remove hypotheses, which violates via isotope annotation discovered ion charge
  if(length(isotopes) > 0){
    hypothese <- checkIsotopes(hypothese, isotopes, ipeak);
  }

  if(nrow(hypothese) < 2){
    return(NULL);
  }

  #Test if hypothese grps include mandatory ions
  #Filter Rules #2
  if(length(quasimolion) > 0){
    hypothese <- checkQuasimolion(hypothese, quasimolion);
  }

  if(nrow(hypothese) < 2){
    return(NULL);
  };

  #Entferne Hypothesen, welche gegen OID-Score&Kausalität verstossen!
  hypothese <- checkOidCausality(hypothese, rules[rules.idx, ]);
  if(nrow(hypothese) < 2){
    return(NULL);
  };

  #Prüfe IPS-Score
  hypothese <- checkIps(hypothese)
  if(nrow(hypothese) < 2){
    return(NULL)
  }

  #We have hypotheses and want to add neutral losses
  if("typ" %in% colnames(rules)){
    hypothese <- addFragments(hypothese, rules, mz)

    hypothese <- resolveFragmentConnections(hypothese)
  }
  return(hypothese);
}


createHypothese <- function(ML, rules, devppm, mzabs, naIdx){
  ML.nrow <- nrow(ML);
  ML.vec <- as.vector(ML);
  max.value <- max(round(ML, 0));
  hashmap <- vector(mode="list", length=max.value);

  print("---Setup hash map")
  for(i in 1:length(ML)){
    val <- trunc(ML[i],0);
    if(val>1){
      hashmap[[val]] <- c(hashmap[[val]],i);
    }
  }
  if("ips" %in% colnames(rules)){
    score <- "ips"
  }else{
    score <- "score"
  }


  hypothese <- matrix(NA,ncol=8,nrow=0);
  colnames(hypothese) <- c("massID", "ruleID", "nmol", "charge", "mass", "score", "massgrp", "check");
  massgrp <- 1;

  print("---Check along hashmap")
  max_length <- length(seq(along=hashmap))
  for(i in seq(along=hashmap)){
    print(paste('------ ', i, 'of', max_length))
    if(is.null(hashmap[[i]])){
      next;
    }

    candidates <- ML.vec[hashmap[[i]]];
    candidates.index <- hashmap[[i]];
    if(i != 1 && !is.null(hashmap[[i-1]]) && min(candidates) < i+(2*devppm*i+mzabs)){
      index <- which(ML.vec[hashmap[[i-1]]]> i-(2*devppm*i+mzabs))
      if(length(index)>0) {
        candidates <- c(candidates, ML.vec[hashmap[[i-1]]][index]);
        candidates.index <- c(candidates.index,hashmap[[i-1]][index]);
      }
    }

    if(length(candidates) < 2){
      next;
    }

    tol <- max(2*devppm*mean(candidates, na.rm=TRUE))+ mzabs;
    result <- cutree(hclust(dist(candidates)), h=tol);
    index <- which(table(result) >= 2);
    if(length(index) == 0){
      next;
    }
    #print(paste("----Check candidates", i, "in", length(hashmap)))

    m <- lapply(index, function(x) which(result == x));
    for(ii in 1:length(m)){
      ini.adducts <- candidates.index[m[[ii]]];
      for( iii in 1:length(ini.adducts)){
        adduct <- ini.adducts[iii] %/% ML.nrow +1;
        mass   <- ini.adducts[iii] %% ML.nrow;
        if(mass == 0){
          mass <- ML.nrow;
          adduct <- adduct -1;
        }
        hypothese <- rbind(hypothese, c(naIdx[mass], adduct, rules[adduct, "nmol"], rules[adduct, "charge"], mean(candidates[m[[ii]]]),  rules[adduct,score],massgrp ,1));
      }
      if(length(unique(hypothese[which(hypothese[, "massgrp"] == massgrp), "massID"])) == 1){
        ##only one mass annotated
        hypothese <- hypothese[-(which(hypothese[,"massgrp"]==massgrp)),,drop=FALSE]
      }else{
        massgrp <- massgrp +1;
      }
    }
  }
  return(hypothese);
}

resolveFragmentConnections <- function(hypothese){
  #Order hypothese after mass
  hypothese <- hypothese[order(hypothese[, "mass"], decreasing=TRUE), ]

  for(massgrp in unique(hypothese[, "massgrp"])){
    index <- which(hypothese[, "massgrp"] == massgrp & !is.na(hypothese[, "parent"]))
    if(length(index) > 0) {
      index2 <- which(hypothese[, "massID"] %in% hypothese[index, "massID"] & hypothese[, "massgrp"] != massgrp)
      if(length(index2) > 0){
        massgrp2del <- which(hypothese[, "massgrp"] %in% unique(hypothese[index2, "massgrp"]))
        hypothese <- hypothese[-massgrp2del, ]
      }
    }
  }
  return(hypothese)
}

addFragments <- function(hypothese, rules, mz){
  #check every hypothese grp
  fragments <- rules[which(rules[, "typ"] == "F"), , drop=FALSE]
  hypothese <- cbind(hypothese, NA);
  colnames(hypothese)[ncol(hypothese)] <- "parent"
  if(nrow(fragments) < 1){
    #no fragment exists in rules
    return(hypothese)
  }

  orderMZ <- cbind(order(mz),order(order(mz)))
  sortMZ <- cbind(mz,1:length(mz))
  sortMZ <- sortMZ[order(sortMZ[,1]),]

  for(massgrp in unique(hypothese[, "massgrp"])){
    for(index in which(hypothese[, "ruleID"] %in% unique(fragments[, "parent"]) &
                       hypothese[, "massgrp"] == massgrp)){
      massID <- hypothese[index, "massID"]
      ruleID <- hypothese[index, "ruleID"]
      indexFrag <- which(fragments[, "parent"] == ruleID)

      while(length(massID) > 0){
        result <- fastMatch(sortMZ[1:orderMZ[massID[1],2],1], mz[massID[1]] +
                              fragments[indexFrag, "massdiff"], tol=0.05)
        invisible(sapply(1:orderMZ[massID[1],2], function(x){
          if(!is.null(result[[x]])){
            massID <<- c(massID, orderMZ[x,1]);
            indexFrags <- indexFrag[result[[x]]];
            tmpRes <- cbind(orderMZ[x,1], as.numeric(rownames(fragments)[indexFrags]), fragments[indexFrags, c("nmol", "charge")],
                            hypothese[index, "mass"], fragments[indexFrags, c("score")],
                            massgrp, 1, massID[1], deparse.level=0)
            colnames(tmpRes) <- colnames(hypothese)
            hypothese <<- rbind(hypothese, tmpRes);
          }
        }))
        massID <- massID[-1];
      }
    }
  }
  return(hypothese)
}


getderivativeIons <- function(annoID, annoGrp, rules, npeaks){
  derivativeIons <- vector("list", npeaks);
  charge <- 0;
  #check that we have annotations
  if(nrow(annoID) < 1){
    return(derivativeIons);
  }

  for(i in 1:nrow(annoID)){
    peakid  <-  annoID[i, 1];
    grpid   <-  annoID[i, 2];
    ruleid  <-  annoID[i, 3];
    parent  <-  annoID[i, 4];
    #         if(is.null(derivativeIons[[peakid]])){
    #           #Peak has no annotation
    #           if(charge == 0 | rules[ruleid, "charge"] == charge){
    #             derivativeIons[[peakid]][[1]] <- list(rule_id = ruleid, charge = rules[ruleid, "charge"],
    #                                                   nmol= rules[ruleid, "nmol"], name=paste(rules[ruleid, "name"]),
    #                                                   mass=annoGrp[which(annoGrp[, "id"] == grpid), 2],parent=parent)
    #             }
    #         }else{
    #           #Peak has already annotations
    if(charge == 0 | rules[ruleid, "charge"] == charge){
      mass <- annoGrp[which(annoGrp[, "id"] == grpid), 2];
      if(is.na(parent)){
        name <- paste(rules[ruleid, "name"]);
      } else {
        #look for name
        name <- paste(rules[ruleid, "name"]);
        for(ii in seq(along=derivativeIons[[parent]])){
          if(derivativeIons[[parent]][[ii]]$mass == mass){

            break;
          }
        }
      }
      derivativeIons[[peakid]][[(length(derivativeIons[[peakid]])+1)]] <-
        list(rule_id = ruleid, charge=rules[ruleid, "charge"],
             nmol=rules[ruleid, "nmol"], name=name,
             mass=mass,
             parent=parent)
    }
    #         }
    charge=0;
  }
  return(derivativeIons);
}

checkIps <- function(hypothese){
  for(hyp in 1:nrow(hypothese)){
    if(length(which(hypothese[, "massgrp"] == hypothese[hyp, "massgrp"])) < 2){
      hypothese[hyp, "check"] = 0;
    }
  }
  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), ];
  if(is.null(nrow(hypothese))) {
    hypothese <- matrix(hypothese, byrow=F, ncol=9)
  }
  if(nrow(hypothese) < 1){
    colnames(hypothese)<-c("massID", "ruleID", "nmol", "charge", "mass", "oidscore", "ips","massgrp", "check")
    return(hypothese)
  }
  for(hyp in 1:nrow(hypothese)){
    if(length(id <- which(hypothese[, "massID"] == hypothese[hyp, "massID"] & hypothese[, "check"] != 0)) > 1){
      masses <- hypothese[id, "mass"]
      nmasses <- sapply(masses, function(x) {
        sum(hypothese[which(hypothese[, "mass"] == x), "score"])
      })
      masses <- masses[-which(nmasses == max(nmasses))];
      if(length(masses) > 0){
        hypothese[unlist(sapply(masses, function(x) {which(hypothese[, "mass"]==x)})), "check"]=0;
      }
    }
  }

  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), ,drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)
  return(hypothese)
}

checkOidCausality <- function(hypothese,rules){
  #check every hypothese grp
  for(hyp in unique(hypothese[,"massgrp"])){
    hyp.nmol <- which(hypothese[, "massgrp"] == hyp & hypothese[, "nmol"] > 1)

    for(hyp.nmol.idx in hyp.nmol){
      if(length(indi <- which(hypothese[, "mass"] == hypothese[hyp.nmol.idx, "mass"] &
                              abs(hypothese[, "charge"]) == hypothese[, "nmol"])) > 1){
        if(hyp.nmol.idx %in% indi){
          #check if [M+H] [2M+2H]... annotate the same molecule
          massdiff <- rules[hypothese[indi, "ruleID"], "massdiff"] /
            rules[hypothese[indi, "ruleID"], "charge"]
          if(length(indi_new <- which(duplicated(massdiff))) > 0){
            hypothese[hyp.nmol.idx, "check"] <- 0;
          }
        }
      }
    }
  }

  #     #check nmol
  #     if(hypothese[hyp, "nmol"] > 1){
  #       #nmol > 1;
  #       checkSure <- TRUE;
  #       if(hypothese[hyp, "charge"] == 1){
  #         #nmol > 1 and charge = 1; e.g. [2M+H]+, ensure [M+H] is there
  #         if(length(which(hypothese[, "mass"] == hypothese[hyp, "mass"] & hypothese[, "oidscore"] == hypothese[hyp, "oidscore"])) > 1){
  #           #same oidscore is there, could also be [3M+H]; otherwise could not check
  #           for(prof in (hypothese[hyp, "nmol"] - 1):1){
  #           #check if [M+H] is there, for a [3M+H], [2M+H] and [M+H] has to be there
  #             indi <- which(hypothese[,"mass"] == hypothese[hyp,"mass"] & hypothese[,"oidscore"] == hypothese[hyp,"oidscore"] & hypothese[,"nmol"] == prof)
  #             if(length(indi) == 0){
  #               checkSure <- FALSE;
  #               hypothese[hyp,"check"] <- 0;
  #               next;
  #             }
  #           }
  #         }
  #       }else if(abs(hypothese[hyp, "charge"]) == hypothese[hyp, "nmol"]){
  #         #nmol > 1 and charge = nmol; e.g. [2M+2H]2+
  #         if(length(which(hypothese[, "mass"] == hypothese[hyp, "mass"] & hypothese[, "oidscore"] == hypothese[hyp, "oidscore"])) > 1){
  #           for(prof in (hypothese[hyp,"nmol"]-1):1){
  #             indi<-which(hypothese[,"mass"]==hypothese[hyp,"mass"] & hypothese[,"oidscore"]== hypothese[hyp,"oidscore"] & hypothese[,"nmol"]==prof)
  #             if(length(indi) == 0){
  #               checkSure <- FALSE;
  #               hypothese[hyp,"check"] <- 0;#next;
  #             }
  #           }
  #         }
  #         if(length(indi <- which(hypothese[, "mass"] == hypothese[hyp, "mass"] & abs(hypothese[, "charge"]) == hypothese[, "nmol"])) > 1){
  #           #check if [M+H] [2M+2H]... annotate the same molecule
  #           massdiff <- rules[hypothese[indi, "ruleID"], "massdiff"] / rules[hypothese[indi, "ruleID"], "charge"]
  #           if(length(indi_new <- which(duplicated(massdiff))) > 0){
  #             checkSure <- FALSE;
  #             hypothese[hyp, "check"] <- 0;
  #           }
  #         }
  #       }
  #       if(checkSure){
  #         hypothese[hyp, "check"] <- 1;
  #       }
  #     }
  #   }
  hypothese <- hypothese[which(hypothese[, "check"] == TRUE), ,drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)
  return(hypothese)
}

checkQuasimolion <- function(hypothese, quasimolion){
  hypomass <- unique(hypothese[, "mass"])
  for(mass in 1:length(hypomass)){
    if(!any(quasimolion %in% hypothese[which(hypothese[, "mass"] == hypomass[mass]), "ruleID"])){
      hypothese[which(hypothese[, "mass"] == hypomass[mass]), "check"] = 0;
    }else if(is.null(nrow(hypothese[which(hypothese[, "mass"] == hypomass[mass]), ]))){
      hypothese[which(hypothese[, "mass"] == hypomass[mass]), "check"] = 0;
    }
  }

  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), , drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)

  return(hypothese)
}

checkIsotopes <- function(hypothese, isotopes, ipeak){
  for(hyp in 1:nrow(hypothese)){
    peakid <- ipeak[hypothese[hyp, 1]];
    if(!is.null(isotopes[[peakid]])){
      #Isotope da
      explainable <- FALSE;
      if(isotopes[[peakid]]$charge == abs(hypothese[hyp, "charge"])){
        explainable <- TRUE;
      }
      if(!explainable){
        #delete Rule
        hypothese[hyp,"check"]=0;
      }
    }
  }
  hypothese <- hypothese[which(hypothese[, "check"]==TRUE), ,drop=FALSE];
  #check if hypothese grps annotate at least two different peaks
  hypothese <- checkHypothese(hypothese)

  return(hypothese)
}

checkHypothese <- function(hypothese){
  if(is.null(nrow(hypothese))){
    hypothese <- matrix(hypothese, byrow=F, ncol=8)
  }
  colnames(hypothese) <- c("massID", "ruleID", "nmol", "charge", "mass", "score", "massgrp", "check")
  for(i in unique(hypothese[,"massgrp"])){
    if(length(unique(hypothese[which(hypothese[, "massgrp"] == i), "massID"])) == 1){
      ##only one mass annotated
      hypothese <- hypothese[-(which(hypothese[,"massgrp"]==i)), , drop=FALSE]
    }
  }
  return(hypothese)
}


add_same_oidscore <-function(hypo,adducts,adducts_no_oid){
  hypo_new<-matrix(NA,ncol=6)
  colnames(hypo_new)<-c("massID","ruleID","nmol","charge","mass","oidscore")
  ids<-hypo[,"ruleID"];
  for(i in 1:nrow(hypo))
  {
    index<-which(adducts[,"oidscore"]==adducts_no_oid[ids[i],"oidscore"])
    hypo_new<-rbind(hypo_new,matrix(cbind(hypo[i,"massID"],index,adducts[index,"nmol"],adducts[index,"charge"],hypo[i,"mass"],adducts[index,"oidscore"]),ncol=6))
  }
  hypo_new<-hypo_new[-1,];
  return(hypo_new);
}


massDiffMatrix <- function(m, rules){
  #m - m/z vector
  #rules - annotation rules
  nRules <- nrow(rules);
  DM   <- matrix(NA, length(m), nRules)

  for (i in seq_along(m)){
    for (j in seq_len(nRules)){
      DM[i, j] <- (abs(rules[j, "charge"] * m[i]) - rules[j, "massdiff"]) / rules[j, "nmol"]    # ((z*m) - add) /n
    }
  }
  return(DM)
}

massDiffMatrixNL <- function(m,neutralloss){
  nadd <- nrow(neutralloss)
  DM <- matrix(NA,length(m),nadd)

  for (i in 1:length(m))
    for (j in 1:nadd)
      DM[i,j] <- m[i] - neutralloss[j,"massdiff"]

  return(DM)
}


naOmit <- function(x) {
  return (x[!is.na(x)]);
}
