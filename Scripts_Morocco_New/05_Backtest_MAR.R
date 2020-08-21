######## SCRIPT MODEL SSE Backtest script #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)




####################################################################################-
#------------------------------------ Backtest  -----------------------------------
####################################################################################-

#Point de départ : 
coefficients_SSE_start <- SSE_coeffcients_males_df[,'2006']

#Point d'arrivée :
coefficients_SSE_end <- SSE_coeffcients_males_df[,'2016']

#Optimisations successives avec 50% hump 50% infant et target = amélioration historique

#Récupération des améliorations historiques
histo_LEimp <- subset(compute_Ex(), sex == 'Male' & year > 2006)
plot_ExImp()
LE_backtest <- c()

coeff_SSE <- coefficients_SSE_start

for (year in histo_LEimp$year){
  optimization_results <- optimize_LE(target_LE = histo_LEimp[as.character(year), 'var'], perc_Infant = 0.5, SSE_coefficients = coeff_SSE)  
  
  coeff_SSE <-compute_newCoeff_optim(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , SSE_coefficients =  coeff_SSE)
  
  LE_backtest[as.character(year)] <- compute_Ex_optim(coeff_SSE)
}


grid.arrange(plot_qx_sensis_optim(SSE_coefficients = coefficients_SSE_start), 
plot_qx_sensis_optim(SSE_coefficients = coeff_SSE), 
plot_qx_sensis_optim(SSE_coefficients = coefficients_SSE_end))



compute_Ex_optim(coeff_SSE) - compute_Ex_optim(coefficients_SSE_end)



#BACKTEST PLOT 

raw_data <- SSE_data_male_df[,'2016']
plot_backtest(coeff_SSE, raw_data)






