library(ggplot2)

full_data = NULL
for(name in list.files(path = "CI_outputs/")[-1]){
  temp_rows = read.csv(paste("CI_outputs/",name, sep =""))
  full_data = rbind(full_data, temp_rows)
}

names(full_data) = c("Lower_CIs", "Upper_CIs", "Iteration", "Knot Value")
full_data$`Knot Value` = as.factor(as.character(full_data$`Knot Value`))
levels(full_data$`Knot Value`) = c("0","0.1750","0.2750","0.3750","0.4750",
                                   "0.5750","0.6750","0.7750","1.0000")

ggplot(full_data, aes(x=Lower_CIs, color=`Knot Value`)) + stat_ecdf() + 
  #scale_color_discrete(name = "description", labels = c("Upper CI on Mu"))+
  xlab("Nominal Coverage")+ylab("Empirical Coverage")+ 
  theme(text = element_text(size=30), axis.text.x = element_text(angle = 90), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))

