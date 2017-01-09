library(foreign)
library(ggplot2)

##############################################
# Load the data set
##############################################
prias = read.spss(file.choose(), to.data.frame=TRUE)
colnames(prias)[match("Discontinued", colnames(prias))] = "DiscontinuedType"
prias$REASO0 = NULL
prias$Dummy.0 = rep(0, nrow(prias))

long_columns = c(c(15:49), c(50:84), c(85:119), c(120:154), c(159:193))

##############################################
# Remove subjects
##############################################
#15 patients do not have any information and they all have Age=NA
empty_patients = prias[is.na(prias$Age),]

#Some patients are too young to be real. So basically remove both empty ones and these patients too
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
prias$P_ID = droplevels(as.factor(prias$P_ID))

# DO NOT ATTEMPT THE COMMENTED. THIS GIVES ISSUE WHEN MERGING THIS COLUMN WITH THE OTHER 35 to make a long dataset
# prias$Gleason_sum = as.ordered(prias$Gleason_sum)
# prias$Gleason1_2 = as.ordered(prias$Gleason1_2)
# prias$Gleason2_2 = as.ordered(prias$Gleason2_2)

# DO NOT DO THIS. COMMENTED NOW. For survival analysis we want it in 0,1 form
# prias$DiscontinuedYesNo = factor(prias$DiscontinuedYesNo, labels = c("No", "Yes"))

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
condition1 = prias$DiscontinuedYesNo==0 & 
  (!is.na(prias$DiscontinuedType) | !is.na(prias$Reason_treatment))

View(prias[condition1, -long_columns])

#Temporary treatment: Remove them
prias = prias[!condition1,]

###########
#Patients who are reported YES discontinued but don't have date of discontinuation
condition2 = prias$DiscontinuedYesNo==1 & is.na(prias$Date_discontinued)
View(prias[condition2, -long_columns])

#Temporary treatment: Remove them
prias = prias[!condition2,]

###########
#Patients who discontinued YES but don't have DiscontinuedType yet they do have reason of treatment
condition3 = prias$DiscontinuedYesNo==1 & 
  is.na(prias$DiscontinuedType) & !is.na(prias$Reason_treatment)
View(prias[condition3, -long_columns])

#Temporary treatment: Remove them because Other, Based on protocol advice could be Lost to follow up or Treatment...we dunno which is reality
prias=prias[!condition3,]

###########
#Patients who discontinued but neither have DiscontinuedType nor have reason of treatment
condition4 = prias$DiscontinuedYesNo==1 & 
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
#Patients who Died, took treatment, watchful waiting BUT reason of treament was lost to follow up
condition6 = prias$Reason_treatment %in% "Lost to FU" & !(prias$DiscontinuedType %in% "Lost to FU")
View(prias[condition6, -long_columns])

#Temporary treament: Do nothing, because
#This is possible if first the reason of discontinuation column was filled and later the discontinuation type
NULL

prias$P_ID = droplevels(prias$P_ID)

#####################################################
# Create simplified versions of various columns
#####################################################
prias$diDRE = ordered(ifelse(prias$DRE %in% c("T1b", "T1c"), yes = "T1", no = "T2"))
prias$diGleason_sum = ordered(ifelse(prias$Gleason_sum<=6, yes="Low", no="High"), levels=c("Low", "High"))

#####################################################
# Categories for discontinuation
#####################################################
#By now, everyone who discontinued has a reason
any(is.na(prias$DiscontinuedYesNo))
any(prias$DiscontinuedYesNo==1 & is.na(prias$DiscontinuedType))

prias$event_type = factor(sapply(1:nrow(prias), function(index){
  prias_i = prias[index, ]
  
  if(prias_i$DiscontinuedYesNo==1){
    if(prias_i$DiscontinuedType %in% "Lost to FU"){
      #Reasons are NA, Other, based on protocol and Lost to FU
      return("Censored")
    }else if(prias_i$DiscontinuedType %in% "Watchful waiting"){
      #Reasons are NA, Anxiety, Based on protocol, lost to FU and other 
      return("Watchful waiting")
    }else if(prias_i$DiscontinuedType %in% "Died"){
      if(is.na(prias_i$Reason_treatment)){
        return("Died-Progression")
      }else{
        #Reasons are Other and Lost to FU
        return("Died-Other")
      }
    }else if(prias_i$DiscontinuedType %in% "On request" | prias_i$Reason_treatment %in% c("Anxiety", "On request", "Other")){
      return("Anxiety/Other/Request")
    }else{
      return("Treatment-Progression")
    }
  }else{
    return(NA)
  }
}, simplify = T))


#################################################
#################################################
# Convert wide to long and order by patient id and time
#################################################
#################################################
prias_long=reshape(prias, direction='long', idvar='P_ID', timevar = "visit_number",
        varying=list(c(3,50:84), c(4,15:49), c(5, 120:154), c(8, 85:119), c(194, 159:193)),
        v.names=c('dom', 'psa', 'dre', 'gleason', 'dummy'))
prias_long = prias_long[order(prias_long$P_ID, prias_long$dom, na.last = T), ]
prias_long$dummy = factor(prias_long$dummy, labels = c("No", "Yes"))

#Why divide by 10 as well. Because there was an extra 10 in the time value. 
prias_long$visitTimeDays = unlist(tapply(prias_long$dom, prias_long$P_ID, function(x){(x-x[1])/(24*60*60*10)}))
prias_long$visitTimeYears = prias_long$visitTimeDays/365

#Some measurements were dummy. i.e. patient didnt turn up. we return psa and dom from sapply
condition_long_1 = prias_long$dummy %in% "Yes"
prias_long$psa[condition_long_1] = NA
prias_long$dre[condition_long_1] = NA
prias_long$gleason[condition_long_1] = NA
prias_long$dom[condition_long_1] = NA

#Gleason scores which are 0 or 1(there is one guy with 1) should be counted as NA
condition_long_2 = prias_long$gleason %in% c(0,1)
prias_long$gleason[condition_long_2] = NA

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

#Cases where psa/gleason/dre is available but date of measurement is missing
condition_long_3 = is.na(prias_long$dom) & (!is.na(prias_long$psa) | !is.na(prias_long$dre) |
                                              !is.na(prias_long$gleason))
View(prias_long[condition_long_3,])
#Temporary solution: remove them for now, otherwise take the date it should have actually been measured
prias_long = prias_long[!condition_long_3,]

#Cases where first two measurements have same DOM but different PSA scores. There are two such people
p_id_dom1dom2same = unique(prias_long$P_ID)[unlist(tapply(prias_long$visitTimeDays, prias_long$P_ID, function(x){sum(x==0, na.rm = T)>1}))]
View(prias_long[prias_long$P_ID %in% p_id_dom1dom2same,])

# No idea how to deal with this one

#Check now if visittimes are sorted for every person. yes they are sorted
any(unlist(tapply(prias_long$visitTimeDays, prias_long$P_ID, function(x){is.unsorted(x, na.rm = T)})))

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

