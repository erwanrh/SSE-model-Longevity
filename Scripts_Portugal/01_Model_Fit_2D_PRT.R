######## SCRIPT MODEL SSE Fit Data and Create dataframes #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)


####################################################################################-
#---------------------------MODEL FIT 2dimensions   ----------------------------
####################################################################################-
## FIT DU MODELE SUR PLUSIEURS ANNEES 
country_code <- 'PRT'
year_min <- 1960
year_max <- 2018

SSE_coeffcients_males <-list()
SSE_coeffcients_females <- list()
SSE_deathrates_male <- list()
SSE_deathrates_female <- list()
SSE_splines_male <- NULL
SSE_splines_female <-NULL
SSE_data_male <- list()
SSE_data_female <- list()

Data_males2D_e <- c()
Data_males2D_x <- c()
Data_males2D_d <- c()

#FIT du modèle sur chaque année entre 1900 et 2017
for (year in seq(2015, 2017)){
  Data_males <- HMD2MH(country=country_code,year=year, sex='males',path='Data',xtra=TRUE)
  Data_males$n[ Data_males$n==0] <-0.01
  Data_males$m[ Data_males$d ==0] <-  0

  
  Data_males2D_e<- cbind(Data_males2D_e, Data_males$n)
  Data_males2D_x <- cbind(Data_males2D_x, Data_males$x)
  Data_males2D_d <- cbind(Data_males2D_d, Data_males$d)
}

Data_males2D_e <- data.frame(Data_males2D_e)
Data_males2D_x <- data.frame(Data_males2D_x)
Data_males2D_d <- data.frame(Data_males2D_d)



fit <- sse2.fit(Data_males2D_x, Data_males2D_d, Data_males2D_e, 100, 1, kappa = 10^6, deg = NULL,
               plotIT = TRUE, plotFIT = TRUE, lab = " ", mon = TRUE,
               max.it=200, ridge=10^-4,
               x1 = 35, x2 = 50)




# DATA FRAMES Coefficients et Death Rates
SSE_coeffcients_males_df <- t(data.frame(matrix(unlist(SSE_coeffcients_males), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_coeffcients_females_df <- t(data.frame(matrix(unlist(SSE_coeffcients_females), nrow=length(SSE_coeffcients_females), byrow=T),row.names = names(SSE_coeffcients_females)))
SSE_deathrates_male_df <-  t(data.frame(matrix(unlist(SSE_deathrates_male), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_deathrates_female_df <-  t(data.frame(matrix(unlist(SSE_deathrates_female), nrow=length(SSE_coeffcients_females), byrow=T),row.names = names(SSE_coeffcients_females)))
SSE_data_male_df <- t(data.frame(matrix(unlist(SSE_data_male), nrow=length(SSE_data_male), byrow=T),row.names = names(SSE_data_male)))
SSE_data_female_df <- t(data.frame(matrix(unlist(SSE_data_female), nrow=length(SSE_data_female), byrow=T),row.names = names(SSE_data_female)))

row.names(SSE_deathrates_male_df) <- 0:110
row.names(SSE_deathrates_female_df) <- 0:110

