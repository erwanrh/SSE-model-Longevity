######## SCRIPT MODEL SSE Fit Data and Create dataframes #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)


####################################################################################-
#---------------------------MODEL FIT MOROCCO ----------------------------
####################################################################################-

Data_qx_MAR <- read.csv('Data_France/MAR/WHO_DATA/Morocco_plot_qx_polyfit.csv', row.names = 1)
Data_dx_MAR <- read.csv('Data_France/MAR/WHO_DATA/Morocco_plot_deaths_polyfit.csv', row.names = 1)

ggplot(Data_dx_MAR) + geom_line(aes(x= Age, y = qx_male, group=year, color=year)) + 
  scale_y_continuous(trans = 'log10')





year_ <- 2016
data_temp <- subset(Data_MAR, sex== 'M' & year == year_, select = c('age', 'deaths', 'expo', 'qx') )
colnames(data_temp) <- c('x','d','n','m')

SSE_males <- morthump(data=data_temp, model='sse')
    
#Récupération des coefficients alpha
SSE_coeffcients_males_df <-as.data.frame(as.matrix(coef(SSE_males)))
    #Récupération des taux fittés
SSE_deathrates_male_df<- as.data.frame(SSE_males$mhat$mhat3[,1]+SSE_males$mhat$mhat1[,1]+SSE_males$mhat$mhat2[,1])
colnames(SSE_deathrates_male_df) <- '2014'

SSE_data_male <- as.data.frame(SSE_males$data$m)





