######## SCRIPT MODEL SSE data for France #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)


####################################################################################-
#--------------------------- Optimisation / Projection ----------------------------
####################################################################################-

compute_delta_LE_months <- function(sensis_infant=0, sensis_hump=0, sensis_senescent=0, period = 2017){

  compute_delta_Ex_sensis(sensis_infant= sensis_infant, sensis_hump= sensis_hump, sensis_senescent= sensis_senescent, period = period) * 12
  
}



optim_LE_Hump <- function(target_LE_Hump, sensis_Hump){
  
  target_LE_Hump - compute_delta_LE_months(sensis_infant = 0, sensis_hump = sensis_Hump, sensis_senescent = 0, period = 2016)
  
}


optim_LE_Infant <- function(target_LE_Infant, sensis_Infant){
  
  target_LE_Infant - compute_delta_LE_months(sensis_infant = sensis_Infant, sensis_hump = 0, sensis_senescent = 0, period = 2016)
  
}



optimize_LE <- function(target_LE, perc_Infant){
  
  #Calcul de la cible d'amélioration pour chaque composante
  target_LE_Infant <- perc_Infant * target_LE
  target_LE_Hump <- (1-perc_Infant) * target_LE
  
  #Solving 
  Infant_root <- uniroot(f = optim_LE_Infant, lower = 0, interval = c(-0,10), target_LE_Infant= target_LE_Infant)$root
  Hump_root <- uniroot(f = optim_LE_Hump, lower = 0, interval = c(-0,10), target_LE_Hump = target_LE_Hump)$root
  
  list('Infant_imp' = Infant_root, 'Hump_imp' = Hump_root)
  
}


