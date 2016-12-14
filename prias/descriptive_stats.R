########################################
# Graphs at baseline
########################################
ggplot(data=prias) + geom_histogram(aes(Age, fill=diG)) + facet_grid(diDRE~DiscontinuedYesNo)
ggplot(data=prias) + geom_histogram(aes(Age, fill=diGleason_sum)) + facet_grid(~DiscontinuedYesNo)

#PSA , Gleason and DRE score spaghetti plot
numSub = 1000
selectedSubPID = sample(prias$P_ID, numSub)
ggplot(data=prias_long[prias_long$P_ID %in% selectedSubPID,], 
       aes(x=time, y=psa, group=P_ID)) + 
  geom_line() +
  stat_summary(aes(group=1), colour="#FD5F00", geom = "point", fun.y = median)+
  ylim(0, 75)+
  facet_grid(~DiscontinuedYesNo)


#Gleason should be plotted against time or biopsy number? I guess biopsy number is a better idea
ggplot(data=prias_long[prias_long$P_ID %in% selectedSubPID & !is.na(prias_long$gleason),]) + 
  geom_line(aes(x=biopsy_number, y=gleason, group=P_ID, color=P_ID)) + guides(color=FALSE) + 
  facet_grid(~DiscontinuedYesNo)

#Gleason and DRE percentage stacked bar charts
colorPalette = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                 "#cab2d6", "#6a3d9a")

ggplot(data = prias_long[!is.na(prias_long$biopsy_number),], aes(factor(biopsy_number))) + 
  geom_bar(aes(fill = gleason), position = "fill") + 
  xlab("Time") + ylab("Percent") + facet_grid(~DiscontinuedYesNo) + 
  scale_fill_manual(values=colorPalette[sample(1:10, 10)])

ggplot(data = prias_long[!is.na(prias_long$biopsy_number),], aes(factor(time))) + 
  geom_bar(aes(fill = gleason), position = "fill") + 
  xlab("Time") + ylab("Percent") + facet_grid(~DiscontinuedYesNo) + 
  scale_fill_manual(values=colorPalette[sample(1:10, 10)])

ggplot(data = prias_long[!is.na(prias_long$dre),], aes(factor(time))) + 
  geom_bar(aes(fill = dre), position = "fill") + 
  xlab("Time") + ylab("Percent") + facet_grid(~DiscontinuedYesNo) + 
  scale_fill_manual(values=colorPalette[sample(1:10, 10)])



prias_sub_dataset = prias_long[!is.na(prias_long$biopsy_number),]
prias_sub_dataset$gleason_risk=sapply(prias_sub_dataset$gleason, function(gleason){
  if(gleason<=6){
    "Low Risk"
  }else if(gleason==7){
    "Medium Risk"
  }else{
    "High Risk"
  }
}, simplify = T)




################
#What age group do we have for our patients
#What should we do with the baby patients? Imputation?
################
ggplot(data=prias, aes(x=Age)) + geom_histogram()

###########
#PSA vs DRE vs Gleason
#most are T1c and Gleason score 6
###########
p=ggplot(data = prias, aes(x=PSA)) + geom_histogram()
p+facet_grid(Gleason_sum~DRE)
#p+facet_wrap(~Gleason_sum)

###########
#Spaghetti plot
#A few patients have very high PSA score
###########
ggplot(data = prias_long[!is.na(prias_long$dummy) & prias_long$dummy==0 & prias_long$psa<100,], aes(x=time, y=psa, group=P_ID)) + 
  geom_line()

##########
#Biopsies per person, impact of less number of biopsies on modeling
##########
table(prias$Gleason_sum)

table(prias_sub_dataset$biopsy_number, prias_sub_dataset$gleason_risk)

ggplot(data = prias_sub_dataset, aes(factor(biopsy_number))) + geom_bar(aes(fill = gleason_risk)) + 
  geom_text(stat='count',aes(label=..count..),vjust=-0.5) + xlab("Biopsy Number") + ylab("Count")

ggplot(data = prias_sub_dataset, aes(factor(biopsy_number))) + geom_bar(aes(fill = gleason_risk), position = "fill") + 
  xlab("Biopsy Number") + ylab("Percent")

#############
#Biopsy spaghetti
#############
ggplot(data = prias_sub_dataset, aes(x=time, y=gleason, group=P_ID)) + 
  geom_line()

##########
#number of visits. The order of time is incorrect. Negative values
#########
ggplot(data=prias, aes(x=Nr_FUvisits)) + geom_histogram() + facet_grid(~DiscontinuedYesNo)

ggplot(data=prias_long[prias_long$time<25,], aes(x=diff_dom)) + geom_histogram() + facet_wrap(~time)
