library(foreign)
library(ggplot2)

prias = read.spss(file.choose(), to.data.frame=TRUE)

prias_long=reshape(prias, direction='long', idvar='P_ID',
        varying=list(c(15:49), c(50:84), c(85:119), c(120:154)),
        v.names=c('psa', 'dom', 'gleason', 'dre'))

prias_long = prias_long[order(prias_long$P_ID), ]


# data set with a positive biopsy number
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

##########
#Results
##########
table(prias_sub_dataset$biopsy_number, prias_sub_dataset$gleason_risk)

ggplot(data = prias_sub_dataset, aes(factor(biopsy_number))) + geom_bar(aes(fill = gleason_risk)) + 
  geom_text(stat='count',aes(label=..count..),vjust=-0.5) + xlab("Biopsy Number") + ylab("Count")

ggplot(data = prias_sub_dataset, aes(factor(biopsy_number))) + geom_bar(aes(fill = gleason_risk), position = "fill") + 
  xlab("Biopsy Number") + ylab("Percent")
