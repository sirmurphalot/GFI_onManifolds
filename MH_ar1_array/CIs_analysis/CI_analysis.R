library(ggplot2)

full_data = data.frame(value = NA, description = NA)
for(name in list.files(path = "../LowerRhoCIs/")[-1]){
  temp_data = read.csv(paste("../LowerRhoCIs/",name, sep =""))
  temp_row = data.frame(value = temp_data[1,1], description = "Lower Rho CI")
  full_data = rbind(full_data, temp_row)
}

for(name in list.files(path = "../UpperRhoCIs/")[-1]){
  temp_data = read.csv(paste("../UpperRhoCIs/",name, sep =""))
  temp_row = data.frame(value = temp_data[1,1], description = "Upper Rho CI")
  full_data = rbind(full_data, temp_row)
}

for(name in list.files(path = "../LowerSdCIs/")[-1]){
  temp_data = read.csv(paste("../LowerSdCIs/",name, sep =""))
  temp_row = data.frame(value = temp_data[1,1], description = "Lower SD CI")
  full_data = rbind(full_data, temp_row)
}

for(name in list.files(path = "../UpperSdCIs/")[-1]){
  temp_data = read.csv(paste("../UpperSdCIs/",name, sep =""))
  temp_row = data.frame(value = temp_data[1,1], description = "Upper SD CI")
  full_data = rbind(full_data, temp_row)
}
full_data = full_data[-1,]

ggplot(full_data, aes(x=value,color=description)) + stat_ecdf() + 
  #scale_color_discrete(name = "description", labels = c("Upper CI on Mu"))+
  xlab("Nominal Coverage")+ylab("Empirical Coverage")+ 
  theme(text = element_text(size=30), axis.text.x = element_text(angle = 90), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17))

summary(full_data[which(full_data$description == "Lower Rho CI"),])
summary(full_data[which(full_data$description == "Upper Rho CI"),])
summary(full_data[which(full_data$description == "Lower SD CI"),])
summary(full_data[which(full_data$description == "Upper SD CI"),])

