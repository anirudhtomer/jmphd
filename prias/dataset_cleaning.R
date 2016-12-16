library(foreign)
library(ggplot2)

##############################################
# Load the data set
##############################################
prias = read.spss(file.choose(), to.data.frame=TRUE)
colnames(prias)[match("Discontinued", colnames(prias))] = "DiscontinuedType"
prias$REASO0 = NULL

long_columns = c(c(15:49), c(50:84), c(85:119), c(120:154), c(159:193))

##############################################
# Remove subjects
##############################################
#15 patients do not have any information and they all have Age=NA
empty_patients = prias[is.na(prias$Age),]

#Some patients are too young to be real
prias = prias[!is.na(prias$Age) & prias$Age > 5,]

##############################################
# Repeat biopsy score update
##############################################
#Some patients had a repeat biopsy, so consider their repeat scores as the correct ones
for(i in 1:nrow(prias)){
  gleason_repeat = prias$Gleason1_2[i] + prias$Gleason2_2[i]
  if(!is.na(gleason_repeat) & gleason_repeat > 0){
    prias$Gleason_sum[i] = gleason_repeat
    prias$Num_cores[i] = prias$Num_cores2[i]
    prias$Num_Cores_PC[i] = prias$Num_Cores_PC2[i]
  }
}

##############################################
# Data type cleaning before conversion to long
##############################################
prias$P_ID = as.factor(prias$P_ID)

prias$Gleason_sum = as.ordered(prias$Gleason_sum)
prias$Gleason1_2 = as.ordered(prias$Gleason1_2)
prias$Gleason2_2 = as.ordered(prias$Gleason2_2)

prias$DiscontinuedYesNo = factor(prias$DiscontinuedYesNo, labels = c("No", "Yes"))

levels(prias$DRE) = trimws(levels(prias$DRE))
prias$DRE = as.ordered(droplevels(prias$DRE))

levels(prias$DiscontinuedType) = trimws(levels(prias$DiscontinuedType))
prias$DiscontinuedType[prias$DiscontinuedType==""] = NA
prias$DiscontinuedType = droplevels(prias$DiscontinuedType)

levels(prias$Reason_treatment) = trimws(levels(prias$Reason_treatment))
prias$Reason_treatment[prias$Reason_treatment %in% c("N/A", "")] = NA
prias$Reason_treatment = droplevels(prias$Reason_treatment)


#################################################
#What to do with these patients (wide format data set)???
#################################################

###########
#Patients who are reported NOT discontinued yet have reasons for discontinuation. No date of discontinuation available
condition1 = prias$DiscontinuedYesNo=="No" & 
  (!is.na(prias$DiscontinuedType) | !is.na(prias$Reason_treatment))

View(prias[condition1, -long_columns])

#Temporary treatment: Remove them
prias = prias[!condition1,]

###########
#Patients who are reported YES discontinued but don't have date of discontinuation
condition2 = prias$DiscontinuedYesNo=="Yes" & is.na(prias$Date_discontinued)
View(prias[condition2, -long_columns])

#Temporary treatment: Remove them
prias = prias[!condition2,]

###########
#Patients who discontinued but don't have DiscontinuedType and do have reason of treatment
condition3 = prias$DiscontinuedYesNo=="Yes" & 
  is.na(prias$DiscontinuedType) & !is.na(prias$Reason_treatment)
View(prias[condition3, -long_columns])

#Temporary treatment: Remove them because Other, Based on protocol advice could be Lost to follow up or Treatment...we dunno which is reality
prias=prias[!condition3,]

###########
#Patients who discontinued but neither have DiscontinuedType nor have reason of treatment
condition4 = prias$DiscontinuedYesNo=="Yes" & 
  is.na(prias$DiscontinuedType) & is.na(prias$Reason_treatment)
View(prias[condition4, -long_columns])

#Temporary treatment: Remove them
prias=prias[!condition4,]

###########
#Patients who have Nr_FUvisits = NA. All of these have no scores for psa, dre or gleason on any of the follow up
condition5 = is.na(prias$Nr_FUvisits)
View(prias[condition5, -long_columns])

#Temporary treatment: Remove them
prias=prias[!condition5,]

###########
#Patients who Died, took treatment, watchful waiting and reason of treament was lost to follow up
#This is possible if first the reason of discontinuation column was filled and later the discontinuationt type

#Temporary treament: Do nothing
NULL

#####################################################
# Create simplified versions of various columns
#####################################################
prias$diDRE = ordered(ifelse(prias$DRE %in% c("T1b", "T1c"), yes = "T1", no = "T2"))
prias$diGleason_sum = ordered(ifelse(prias$Gleason_sum<=6, yes="Low", no="High"), levels=c("Low", "High"))

#by now everyone has a discontinued type
prias$event_type = rep("Treatment", nrow(prias))
prias$event_type[prias$DiscontinuedType %in% c("Lost to FU")] = "Censored"
prias$event_type[prias$DiscontinuedType %in% c("Died")] = "Death"

################################################
#Convert wide to long and order by patient id and time
################################################
prias_long=reshape(prias, direction='long', idvar='P_ID', timevar = "visit_number",
        varying=list(c(15:49), c(50:84), c(85:119), c(120:154), c(159:193)),
        v.names=c('psa', 'dom', 'gleason', 'dre', 'dummy'))
prias_long = prias_long[order(prias_long$P_ID, prias_long$dom, na.last = T), ]
prias_long$visit_number = rep(1:35, length(prias$P_ID))
prias_long$dummy = factor(prias_long$dummy, labels = c("No", "Yes"))

#Some measurements were dummy. i.e. patient didnt turn up. we return psa and dom from sapply
prias_long[, c("psa", "dom")] = t(sapply(1:nrow(prias_long), FUN = function(rowNum){
  temp = prias_long[rowNum,]
  if(!is.na(temp$dummy) & temp$dummy=="No"){
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

#########################################################
# Data type cleaning for the long version of the data set
#########################################################
levels(prias_long$dre) = trimws(levels(prias_long$dre))
prias_long$dre = as.ordered(droplevels(prias_long$dre))

#####################################################
# Create simplified versions of various columns
#####################################################
prias_long$didre = rep(NA, nrow(prias_long))
prias_long[!is.na(prias_long$dre), ]$didre = ifelse(prias_long[!is.na(prias_long$dre), ]$dre %in% c("T1c"), yes = "T1", no = ">T1")
prias_long$didre = ordered(prias_long$didre, levels=c("T1", ">T1"))

prias_long$digleason = rep(NA, nrow(prias_long))
prias_long[!is.na(prias_long$gleason), ]$digleason = ifelse(prias_long[!is.na(prias_long$gleason), ]$gleason<=6, yes="Low", no="High")
prias_long$digleason = ordered(prias_long$digleason, levels=c("Low", "High"))

prias_long$gleason = as.ordered(prias_long$gleason)

#################################################
#What to do with these patients (long format data set)???
#################################################

##########
#Cases where psa/gleason/dre is available but date of measurement is missing

View(prias_long[is.na(prias_long$Date_discontinued) & prias_long$DiscontinuedYesNo %in% c("Yes") &
                  (!is.na(prias_long$psa) | !is.na(prias_long$dre) |
                     !is.na(prias_long$gleason)),])
View(prias_long[is.na(prias_long$dom) & (!is.na(prias_long$psa) | !is.na(prias_long$dre) |
                                           !is.na(prias_long$gleason)),])

#Temporary solution: do nothing...take mean of time using information from other patients


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

