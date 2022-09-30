library(ggplot2)
library(reshape2)

full_data = data.frame(value = NA, description = NA)
full_data = read.csv(paste("../CI_outputs/",list.files(path = "../CI_outputs/")[1], sep =""))
for(name in list.files(path = "../CI_outputs/")[-1]){
  temp_data = read.csv(paste("../CI_outputs/",name, sep =""))
  full_data = rbind(full_data, temp_data)
  # temp_row = data.frame(value = temp_data[1,1], description = "Rho Two-Sided CI")
  # full_data = rbind(full_data, temp_row)
}

# For future people who are confused about my code: my matlab array simulation had lower CB and upper CB mixed up.
# I'm correcting it here so I do not have to re-run the simulation.
names(full_data)[1:7] = c("Upper CB, Knot 2", "Upper CB, Knot 3", "Upper CB, Knot 4", "Upper CB, Knot 5", 
                            "Upper CB, Knot 6", "Upper CB, Knot 7", "Upper CB, Knot 8")
temp_full_data = melt(full_data[,1:7])
names(temp_full_data)
ggplot(temp_full_data, aes(x=value,color=variable)) + stat_ecdf() + 
  #scale_color_discrete(name = "description", labels = c("Upper CI on Mu"))+
  xlab("Nominal Coverage")+ylab("Empirical Coverage")+ theme_classic()+ 
  theme(text = element_text(size=30), axis.text.x = element_text(angle = 90), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=10))+
  theme(text=element_text(size=15,  family="sans"))