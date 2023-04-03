library(ggplot2)

full_data = data.frame(value = NA, description = NA)
for(name in list.files(path = "../fixed_rho_outputs/")){
  temp_data = read.csv(paste("../fixed_rho_outputs/",name, sep =""))
  # temp_row = data.frame(value = temp_data[1,3], description = "Lower Rho CI")
  # full_data = rbind(full_data, temp_row)
  temp_row = data.frame(value = temp_data[1,2], description = "Upper Rho CB")
  full_data = rbind(full_data, temp_row)
  # temp_row = data.frame(value = temp_data[1,1], description = "Rho Two-Sided CI")
  # full_data = rbind(full_data, temp_row)
}

for(name in list.files(path = "../fixed_sd_outputs/")){
  temp_data = read.csv(paste("../fixed_sd_outputs/",name, sep =""))
  # temp_row = data.frame(value = temp_data[1,3], description = "Lower SD CI")
  # full_data = rbind(full_data, temp_row)
  temp_row = data.frame(value = temp_data[1,2], description = "Upper SD CB")
  full_data = rbind(full_data, temp_row)
}

full_data = full_data[-1,]

ggplot(full_data, aes(x=value,linetype=description)) + stat_ecdf() + 
  #scale_color_discrete(name = "description", labels = c("Upper CI on Mu"))+
  xlab("Nominal Coverage")+ylab("Empirical Coverage")+ theme_classic()+ 
  theme(text = element_text(size=30), axis.text.x = element_text(angle = 90), axis.title.x = element_text(size=20), 
        axis.title.y = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=17)) +
  theme(text=element_text(size=15,  family="sans")) + geom_hline(yintercept=0.95) + geom_vline(xintercept=c(0.95,1)) + xlim(0,1) +
  ylim(0,1) 

summary(full_data[which(full_data$description == "Lower Rho CI"),])
summary(full_data[which(full_data$description == "Upper Rho CI"),])
summary(full_data[which(full_data$description == "Rho Two-Sided CI"),])
summary(full_data[which(full_data$description == "Lower SD CI"),])
summary(full_data[which(full_data$description == "Upper SD CI"),])
summary(full_data[which(full_data$description == "Lower Sigma CI"),])

mean(full_data[which(full_data$description == "Upper Sigma CI"),]$value >=0.05)
mean(full_data[which(full_data$description == "Upper SD CI"),]$value >=0.1)
mean(full_data[which(full_data$description == "Upper Rho CI"),]$value >=0.1)

