######## SCRIPT MODEL SSE Fit Data and Create dataframes #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)


####################################################################################-
#---------------------------MODEL FIT  ----------------------------
####################################################################################-
## FIT DU MODELE SUR PLUSIEURS ANNEES 
country_code <- 'FR'
year_min <- 1960
year_max <- 2017
step <- 1

SSE_coeffcients_males <-list()
SSE_coeffcients_females <- list()
SSE_deathrates_male <- list()
SSE_deathrates_female <- list()
SSE_splines_male <- NULL
SSE_splines_female <-NULL
SSE_data_male <- list()
SSE_data_female <- list()


#FIT du modèle sur chaque année entre 1900 et 2017
for (year in seq(year_min, year_max, step)){
  print(paste0('Fitting year : ', year, '...'), quote = FALSE) 
  
  #Fetching data from HMD files
  Data_males <- HMD2MH(country=country_code,year=year, sex='males',path='Data',xtra=TRUE)
  Data_females <- HMD2MH(country=country_code,year=year, sex='females',path='Data',xtra=TRUE)
  
  #Data corrections
  # Expo
  Data_males$n[ Data_males$n==0] <-0.01
  Data_females$n[ Data_females$n==0] <-0.01
  #Deaths
  Data_females$d <- as.integer(Data_females$d)
  Data_males$d <- as.integer(Data_males$d)
  #Rates
  Data_females$m[ Data_females$d ==0] <-  0
  Data_males$m[ Data_males$d ==0] <-  0


  if (!any(is.na(Data_males)) & !any(is.na(Data_females)) ){
    
    #Fit
    try(SSE_males<- morthump(data=Data_males, model='sse', x1 = 25, x2 = 40,lambda.sen = 3), silent = TRUE)
    try(SSE_females<- morthump(data=Data_females, model='sse',  x1 = 25, x2 = 40, lambda.sen = 3), silent = TRUE)
    
    
    #Récupération des coefficients alpha
    SSE_coeffcients_males[as.character(year)] <-as.data.frame(as.matrix(coef(SSE_males)))
    SSE_coeffcients_females[as.character(year)] <- as.data.frame(as.matrix(coef(SSE_females)))
    
    #Récupération des taux fittés
    SSE_deathrates_male[as.character(year)] <- as.data.frame(SSE_males$mhat$mhat3[,1]+SSE_males$mhat$mhat1[,1]+SSE_males$mhat$mhat2[,1])
    SSE_deathrates_female[as.character(year)] <- as.data.frame(SSE_females$mhat$mhat3[,1]+SSE_females$mhat$mhat1[,1]+SSE_females$mhat$mhat2[,1])
    
    #Récupération des composantes
    temp_SSE_splines_male <- melt(cbind(SSE_males$XX$X1,SSE_males$XX$X2,SSE_males$XX$X3), value.name = 'splinevalue',varnames = c('age','splinenb'))
    temp_SSE_splines_male$year <- year
    
    temp_SSE_splines_female <- melt(cbind(SSE_females$XX$X1,SSE_females$XX$X2,SSE_females$XX$X3), value.name = 'splinevalue',varnames = c('age','splinenb'))
    temp_SSE_splines_female$year <- year
    
    SSE_splines_male<- rbind(SSE_splines_male,temp_SSE_splines_male)
    SSE_splines_female<- rbind(SSE_splines_female,temp_SSE_splines_female)
    
    SSE_data_male[as.character(year)] <- as.data.frame(SSE_males$data$m)
    SSE_data_female[as.character(year)] <-as.data.frame( SSE_females$data$m)
  }
  
}

# DATA FRAMES Coefficients et Death Rates
SSE_coeffcients_males_df <- t(data.frame(matrix(unlist(SSE_coeffcients_males), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_coeffcients_females_df <- t(data.frame(matrix(unlist(SSE_coeffcients_females), nrow=length(SSE_coeffcients_females), byrow=T),row.names = names(SSE_coeffcients_females)))
SSE_deathrates_male_df <-  t(data.frame(matrix(unlist(SSE_deathrates_male), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_deathrates_female_df <-  t(data.frame(matrix(unlist(SSE_deathrates_female), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_data_male_df <- t(data.frame(matrix(unlist(SSE_data_male), nrow=length(SSE_data_male), byrow=T),row.names = names(SSE_data_male)))
SSE_data_female_df <- t(data.frame(matrix(unlist(SSE_data_female), nrow=length(SSE_data_female), byrow=T),row.names = names(SSE_data_female)))

row.names(SSE_deathrates_male_df) <- 0:110
row.names(SSE_deathrates_female_df) <- 0:110

source('study.R')


