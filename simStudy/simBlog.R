t = 1:10

trajectory = c(1, 1, 1, 4, 4, 4, 4, 10, 10, 10)

ggplot(data = data.frame(Time=t, Biomarker_LVCF=trajectory, Biomarker_LM=fitted(lm(trajectory~t + I(t^2)+ I(t^3)+ I(t^4))))) + 
  geom_step(aes(x=Time, y=Biomarker_LVCF, colour="LVCF"), size=1) + scale_x_continuous(breaks = 1:10) + 
  geom_smooth(aes(x=Time, y=Biomarker_LM, colour="Two stage"), se = F, size=1) + 
  ylab("Biomarker")  +
  scale_colour_manual(name="Method", values=c(LVCF="red", "Two stage"="darkgreen"))


ggplot(data = data.frame(Time=t, Biomarker_LVCF=trajectory, Biomarker_LM=fitted(lm(trajectory~t + I(t^2)+ I(t^3)+ I(t^4))))) + 
  geom_step(aes(x=Time, y=Biomarker_LVCF, colour="LVCF"), size=1) + scale_x_continuous(breaks = 1:10) + 
  geom_smooth(aes(x=Time, y=Biomarker_LM, colour="Joint model"), se = T, size=1) + 
  ylab("Biomarker")  +
  scale_colour_manual(name="Method", values=c(LVCF="red", "Joint model"="blue"))