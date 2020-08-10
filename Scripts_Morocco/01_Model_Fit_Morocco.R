######## SCRIPT MODEL SSE Fit Data and Create dataframes #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)



####################################################################################-
#------------------------------------  FIT MOROCCO ----------------------------------
####################################################################################-

All_interpolated_data <- merge(x = merge(x = Data_nx_MAR_All, y = Data_qx_MAR_All, by = c('Age','year','sex'), suffixes = c('nx', 'qx')),
                               y = Data_dx_MAR_All, by = c('Age','year','sex'), suffixes = c('x', 'dx'))




SSE_coeffcients_males <-list()
SSE_coeffcients_females <- list()
SSE_deathrates_male <- list()
SSE_deathrates_female <- list()
SSE_splines_male <- NULL
SSE_splines_female <-NULL
SSE_data_male <- list()
SSE_data_female <- list()


for (year_ in unique(Data_qx_MAR_All$year)){
  print(paste0('Fitting year ', year_))
  Data_males <- subset(All_interpolated_data, sex== 'M' & year == year_, select = c('Age', 'qx_male', 'qx_malenx', 'qx_maleqx') )
  colnames(Data_males) <- c('x','d','n','m')
  Data_males <- Data_males[order(Data_males$x),]
  
  Data_females <- subset(All_interpolated_data, sex== 'F' & year == year_, select = c('Age', 'qx_male', 'qx_malenx', 'qx_maleqx') )
  colnames(Data_females) <- c('x','d','n','m')
  Data_females <- Data_females[order(Data_females$x),]

  SSE_males<- morthump(data=Data_males, model='sse', lambda.hump = 1, lambda.sen = 1, x1 = 25)
  SSE_females<- morthump(data=Data_females, model='sse', lambda.hump = 1, lambda.sen = 1, x1 = 25)
  
  
  #Récupération des coefficients alpha
  SSE_coeffcients_males[as.character(year_)] <-as.data.frame(as.matrix(coef(SSE_males)))
  SSE_coeffcients_females[as.character(year_)] <- as.data.frame(as.matrix(coef(SSE_females)))
  
  #Récupération des taux fittés
  SSE_deathrates_male[as.character(year_)] <- as.data.frame(SSE_males$mhat$mhat3[,1]+SSE_males$mhat$mhat1[,1]+SSE_males$mhat$mhat2[,1])
  SSE_deathrates_female[as.character(year_)] <- as.data.frame(SSE_females$mhat$mhat3[,1]+SSE_females$mhat$mhat1[,1]+SSE_females$mhat$mhat2[,1])
  

  SSE_data_male[as.character(year_)] <- as.data.frame(SSE_males$data$m)
  SSE_data_female[as.character(year_)] <-as.data.frame( SSE_females$data$m)
  
}






SSE_coeffcients_males_df <- t(data.frame(matrix(unlist(SSE_coeffcients_males), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_coeffcients_females_df <- t(data.frame(matrix(unlist(SSE_coeffcients_females), nrow=length(SSE_coeffcients_females), byrow=T),row.names = names(SSE_coeffcients_females)))
SSE_deathrates_male_df <-  t(data.frame(matrix(unlist(SSE_deathrates_male), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_deathrates_female_df <-  t(data.frame(matrix(unlist(SSE_deathrates_female), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_data_male_df <- t(data.frame(matrix(unlist(SSE_data_male), nrow=length(SSE_data_male), byrow=T),row.names = names(SSE_data_male)))
SSE_data_female_df <- t(data.frame(matrix(unlist(SSE_data_female), nrow=length(SSE_data_female), byrow=T),row.names = names(SSE_data_female)))

row.names(SSE_deathrates_male_df) <- 0:85
row.names(SSE_deathrates_female_df) <- 0:85

