library(foreign)
library(ggplot2)

#Load the data set
prias = read.spss(file.choose(), to.data.frame=TRUE)

colnames(prias)[match("Discontinued", colnames(prias))] = "DiscontinuedType"
prias$REASO0 = NULL
prias$DiscontinuedYesNo = factor(prias$DiscontinuedYesNo, labels = c("No", "Yes"))

#15 patients do not have any information.
#DRE column has blank values and not NA, so it can't be tested
#P_ID column always has values and DiscontinuedYesNo also has values
#Its a bit messy so I just test for Age. Gives same result anyway
empty_patients = prias[is.na(prias$Age),]
prias = prias[!is.na(prias$Age) & prias$Age > 5,]

#Some patients had a repeat biopsy, so consider their repeat scores as the correct ones
prias$gleason_repeat = prias$Gleason1_2 + prias$Gleason2_2
for(i in 1:nrow(prias)){
  if(prias$gleason_repeat[i] > 0){
    prias$Gleason_sum[i] = prias$gleason_repeat[i]
    prias$Num_cores[i] = prias$Num_cores2[i]
    prias$Num_Cores_PC[i] = prias$Num_Cores_PC2[i]
  }
}

#Change data types. Most of them we will change after reshaping to long
prias$P_ID = as.factor(prias$P_ID)
prias$Gleason_sum = as.ordered(prias$Gleason_sum)
prias$Gleason1_2 = as.ordered(prias$Gleason1_2)
prias$Gleason2_2 = as.ordered(prias$Gleason2_2)
prias$gleason_repeat = as.ordered(prias$gleason_repeat)

#Convert wide to long and order by patient id
prias_long=reshape(prias, direction='long', idvar='P_ID',
        varying=list(c(15:49), c(50:84), c(85:119), c(120:154), c(159:193)),
        v.names=c('psa', 'dom', 'gleason', 'dre', 'dummy'))
prias_long = prias_long[order(prias_long$P_ID, prias_long$dom, na.last = T), ]
prias_long$time = rep(1:35, length(prias$P_ID))

#Drop unused levels of DRE and trim the level names
levels(prias$DRE) = trimws(levels(prias$DRE))
levels(prias_long$DRE) = trimws(levels(prias_long$DRE))
prias$DRE = as.ordered(droplevels(prias$DRE))
prias_long$DRE = as.ordered(droplevels(prias_long$DRE))

levels(prias_long$dre) = trimws(levels(prias_long$dre))
prias_long$dre = as.ordered(droplevels(prias_long$dre))

#Rename levels for dummy
prias_long$dummy = factor(prias_long$dummy, labels = c("No", "Yes"))

#Rename levels for DiscontinuedType
levels(prias$DiscontinuedType) = trimws(levels(prias$DiscontinuedType))
levels(prias_long$DiscontinuedType) = trimws(levels(prias_long$DiscontinuedType))

prias_long$DiscontinuedType[prias_long$DiscontinuedType==""] = NA
prias$DiscontinuedType[prias$DiscontinuedType==""] = NA

prias_long$DiscontinuedType = droplevels(prias_long$DiscontinuedType)
prias$DiscontinuedType = droplevels(prias$DiscontinuedType)

#create new type of discontinuation reasons
prias_long$event_type = rep("Treatment", nrow(prias_long))
prias_long$event_type[prias$DiscontinuedType %in% c("Watchful waiting", "Lost to FU", "On request")] = "Censored"
prias_long$event_type[prias$DiscontinuedType %in% c("Watchful waiting", "Lost to FU")] = "Censored"

#Rename levels for Reason treatment
levels(prias$Reason_treatment) = trimws(levels(prias$Reason_treatment))
levels(prias_long$Reason_treatment) = trimws(levels(prias_long$Reason_treatment))

prias_long$Reason_treatment[prias_long$Reason_treatment %in% c("N/A", "")] = NA
prias$Reason_treatment[prias$Reason_treatment %in% c("N/A", "")] = NA

prias_long$Reason_treatment = droplevels(prias_long$Reason_treatment)
prias$Reason_treatment = droplevels(prias$Reason_treatment)

#Some measurements were dummy. i.e. patient didnt turn up. we return psa and dom from sapply
prias_long[, c("psa", "dom")] = t(sapply(1:nrow(prias_long), FUN = function(rowNum){
  temp = prias_long[rowNum,]
  if(!is.na(temp$dummy) && temp$dummy=="No"){
    return(c(temp$psa, temp$dom))
  }else{
    return(c(NA, NA))
  }
}))

#Gleason scores which are 0 should be counted as NA
prias_long$gleason = unlist(lapply(prias_long$gleason, FUN = function(gleason){
  if(!is.na(gleason) && gleason==0){
    return(NA)
  }else{
    return(gleason)
  }
}))

#Baseline gleason scores which are 0 should be counted as NA
prias_long$Gleason_sum = unlist(lapply(prias_long$Gleason_sum, FUN = function(Gleason_sum){
  if(!is.na(Gleason_sum) && Gleason_sum==0){
    return(NA)
  }else{
    return(Gleason_sum)
  }
}))

#PSA scores which are 0 should be counted as NA
prias_long$psa = unlist(lapply(prias_long$psa, FUN = function(psa){
  if(!is.na(psa) && psa==0){
    return(NA)
  }else{
    return(psa)
  }
}))

#How balanced is the data set. as in how often measurements are taken
prias_long$diff_dom = unlist(lapply(prias$P_ID, FUN=function(id){
  tempds = prias_long[prias_long$P_ID==id,]
  numvisits = tempds$Nr_FUvisits[1]
  if(is.na(numvisits)){ #699 such people exist
    rep(NA, 35)
  }else{
    c(0,diff(tempds$dom[1:numvisits]), rep(NA, 35-numvisits))/(24*60*60*10)
    #the extra 10 because there was an extra zero
  }
}))


# Every biopsy for a patient is given a biopsy number
prias_long$biopsy_number = unlist(lapply(prias$P_ID, FUN = function(id){
  biopsy_indicator= as.numeric(prias_long[prias_long$P_ID==id, ]$gleason>0)
  count = 1
  for(i in 1:35){
    if(!is.na(biopsy_indicator[i]) & biopsy_indicator[i]>0){
      biopsy_indicator[i] <- count
      count = count + 1
    }else{
      biopsy_indicator[i] = NA
    }
  }
  biopsy_indicator
}))

#Change data types
prias_long$gleason = as.ordered(prias_long$gleason)

