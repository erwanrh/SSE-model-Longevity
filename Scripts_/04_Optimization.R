######## SCRIPT MODEL SSE Optimization for Morocco #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)


####################################################################################-
#--------------------------- Optimisation / Projection ----------------------------
####################################################################################-
## DELTA EX FOR SENSITIVITY
compute_fitted_sensis_optim <- function(sensis_infant=0, sensis_hump=0, sensis_senescent=0, SSE_coefficients){
  SSE_splines <- read.csv('splinesMAR.csv', row.names = 1)
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(as.matrix(SSE_splines[, 1:2])%*%(as.matrix(SSE_coefficients[1:2]))*(1+sensis_infant)), 
                                               exp(as.matrix(SSE_splines[, 3:27])%*%(SSE_coefficients[3:27])*(1+sensis_senescent)), 
                                               exp(as.matrix(SSE_splines[, 28:52])%*%(SSE_coefficients[28:52])*(1+sensis_hump))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(as.matrix(SSE_splines[, 1:2])%*%(SSE_coefficients[1:2])), 
                                       exp(as.matrix(SSE_splines[, 3:27])%*%(SSE_coefficients[3:27])), 
                                       exp(as.matrix(SSE_splines[, 28:52])%*%(SSE_coefficients[28:52]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <- 0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <- 0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}


compute_delta_Ex_sensis_optim <- function (sensis_infant=0, sensis_hump=0, sensis_senescent=0, SSE_coefficients){
  
  Fitted_sensis_list <- compute_fitted_sensis_optim(sensis_infant, sensis_hump, sensis_senescent, SSE_coefficients)
  Fitted_DataFrame_Sensis1 <- Fitted_sensis_list[['df_sensis']]
  Fitted_DataFrame <- Fitted_sensis_list[['df_base']]
  
  12 * (LE_period_Model(QxModel = Fitted_DataFrame_Sensis1,Age =  0, Period = 'FittedCurve') - LE_period_Model(QxModel = Fitted_DataFrame,Age =  0, Period = 'FittedCurve'))
  
}




optim_LE_Hump <- function(target_LE_Hump, sensis_Hump, SSE_coefficients){
  
  target_LE_Hump - compute_delta_Ex_sensis_optim(sensis_infant = 0, sensis_hump = sensis_Hump, sensis_senescent = 0, SSE_coefficients)
  
}


optim_LE_Infant <- function(target_LE_Infant, sensis_Infant, SSE_coefficients){
  
  target_LE_Infant - compute_delta_Ex_sensis_optim(sensis_infant = sensis_Infant, sensis_hump = 0, sensis_senescent = 0, SSE_coefficients)
  
}



optimize_LE <- function(target_LE, perc_Infant, SSE_coefficients ){
  
  #Calcul de la cible d'amélioration pour chaque composante
  target_LE_Infant <- perc_Infant * target_LE
  target_LE_Hump <- (1-perc_Infant) * target_LE
  
  #Solving 
  Infant_root <- uniroot(f = optim_LE_Infant, lower = 0, interval = c(-0,10), target_LE_Infant= target_LE_Infant, SSE_coefficients = SSE_coefficients)$root
  Hump_root <- uniroot(f = optim_LE_Hump, lower = 0, interval = c(-0,10), target_LE_Hump = target_LE_Hump, SSE_coefficients = SSE_coefficients)$root
  
  list('Infant_imp' = Infant_root, 'Hump_imp' = Hump_root)
  
}


plot_qx_sensis_optim <- function(sensis_infant=0, sensis_hump=0, sensis_senescent=0, SSE_coefficients){
  
  Fitted_sensis_list <- compute_fitted_sensis_optim(sensis_infant, sensis_hump, sensis_senescent, SSE_coefficients)
  Fitted_DataFrame_Sensis1 <- Fitted_sensis_list[['df_sensis']]
  Fitted_DataFrame <- Fitted_sensis_list[['df_base']]
  
  Fitted_DataFrame_plot_sensis1 <- melt(Fitted_DataFrame_Sensis1, id.vars = c('Sensis', 'Age'), variable.name = 'Component', value.name = 'qx')
  Fitted_DataFrame_plot <- melt(Fitted_DataFrame, id.vars = c('Sensis','Age'), variable.name = 'Component', value.name = 'qx')
  
  Final_DF_plot <- rbind(Fitted_DataFrame_plot, Fitted_DataFrame_plot_sensis1)
  
  p<- ggplot(Final_DF_plot) + geom_line(aes(x = Age, y = qx, group = interaction(Sensis, Component), color = Component, linetype = Sensis)) +
    scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ') +
    ggtitle(paste0('Components, Data and Fitted Curve - ', country_code, ' - SSE'))
  
  
  
  
}


