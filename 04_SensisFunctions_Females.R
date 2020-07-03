######## SCRIPT MODEL SSE Sensitivities Analysis #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)




####################################################################################-
#-------------------------- FUNCTIONS SENSITIVITES  ----------------------------
####################################################################################-


compute_fitted_sensis <- function(sensis_infant=0, sensis_hump=0, sensis_senescent=0, period = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period)])*(1+sensis_infant)), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period)])*(1+sensis_senescent)), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period)])*(1+sensis_hump))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <- 0:110
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <- 0:110
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}



#PLOT CURVES BEFORE AND AFTER



## DELTA EX FOR SENSITIVITY
compute_delta_Ex_sensis <- function (sensis_infant=0, sensis_hump=0, sensis_senescent=0, period = 2017){
  
  Fitted_sensis_list <- compute_fitted_sensis(sensis_infant, sensis_hump, sensis_senescent, period)
  Fitted_DataFrame_Sensis1 <- Fitted_sensis_list[['df_sensis']]
  Fitted_DataFrame <- Fitted_sensis_list[['df_base']]
  
  LE_period_Model(QxModel = Fitted_DataFrame_Sensis1,Age =  0, Period = 'FittedCurve') - LE_period_Model(QxModel = Fitted_DataFrame,Age =  0, Period = 'FittedCurve')
  
}

## DELTA entropy FOR SENSITIVITY
compute_delta_entropy_sensis <- function (sensis_infant=0, sensis_hump=0, sensis_senescent=0, period= 2017){
  
  Fitted_sensis_list <- compute_fitted_sensis(sensis_infant, sensis_hump, sensis_senescent, period)
  Fitted_DataFrame_Sensis1 <- Fitted_sensis_list[['df_sensis']]
  Fitted_DataFrame <- Fitted_sensis_list[['df_base']]
  
  Entropy <- subset(compute_entropy_period(Qx.Table = Fitted_DataFrame_Sensis1[,1:4], start.age = 0), period == 'FittedCurve',
                    select = 'value')[1, 'value'] 
  Entropy_Before <- subset(compute_entropy_period(Qx.Table = Fitted_DataFrame[,1:4], start.age = 0), period == 'FittedCurve',
                           select = 'value')[1, 'value'] 
  
  Entropy - Entropy_Before
  
}

plot_qx_sensis <- function (sensis_infant=0, sensis_hump=0, sensis_senescent=0, period = 2017){
  
  Fitted_sensis_list <- compute_fitted_sensis(sensis_infant, sensis_hump, sensis_senescent, period)
  Fitted_DataFrame_Sensis1 <- Fitted_sensis_list[['df_sensis']]
  Fitted_DataFrame <- Fitted_sensis_list[['df_base']]
  
  Fitted_DataFrame_plot_sensis1 <- melt(Fitted_DataFrame_Sensis1, id.vars = c('Sensis', 'Age'), variable.name = 'Component', value.name = 'qx')
  Fitted_DataFrame_plot <- melt(Fitted_DataFrame, id.vars = c('Sensis','Age'), variable.name = 'Component', value.name = 'qx')
  
  Final_DF_plot <- rbind(Fitted_DataFrame_plot, Fitted_DataFrame_plot_sensis1)
  
  p<- ggplot(Final_DF_plot) + geom_line(aes(x = Age, y = qx, group = interaction(Sensis, Component), color = Component, linetype = Sensis)) +
    scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ') +
    ggtitle(paste0('Components, Data and Fitted Curve for France ', period, ' - SSE'))
  
  p
}

## DELTA EX FOR SENSITIVITY
plot_lx_sensis <- function (sensis_infant=0, sensis_hump=0, sensis_senescent=0, period = 2017){
  
  Fitted_sensis_list <- compute_fitted_sensis(sensis_infant, sensis_hump, sensis_senescent, period)
  Fitted_DataFrame_Sensis1 <- Fitted_sensis_list[['df_sensis']]
  Fitted_DataFrame <- Fitted_sensis_list[['df_base']]
  
  
  Fitted_qx <- data.frame(cbind(Fitted_DataFrame$FittedCurve, Fitted_DataFrame_Sensis1$FittedCurve))
  colnames(Fitted_qx) <- c('Base', 'Sensis')
  
  #### Survival function
  lx <- compute_lx_period(Fitted_qx, start.age = 0)
  lx_plot <- melt(lx, varnames = c('Age', 'Period'), value.name = 'lx')
  
  
  
  entropy <- compute_entropy_period(Fitted_qx, start.age = 0)
  entropy$LE <- c(LE_period_Model(0,'Base', QxModel = Fitted_qx),
                  LE_period_Model(0,'Sensis', QxModel = Fitted_qx))
  
  #Graph of fitted survival curves + entropy values
  lx_H_plot <- ggplot(lx_plot) + geom_line(aes(x=Age, y = lx, group = Period, linetype = Period)) +
    annotation_custom(tableGrob(cbind(as.character(entropy[,1]), round(entropy[,2]*100,3),round(entropy[,5],2) ), 
                                rows = NULL, theme=ttheme_default(base_size = 7), cols = c('Period', 'Entropy','LifeEx')),
                      xmin=0, xmax=30, ymin=20000, ymax=25000) +
    ggtitle('Survival curves and entropy on fitted values')
  
  lx_H_plot
}



# Crossed Sensitivities ---------------------------------------------------

compute_fitted_crossedsensis_Senescent <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period2)])), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period1)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period1)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <- 0:110
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <- 0:110
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

compute_fitted_crossedsensis_Hump <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period1)])), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period2)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period1)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <- 0:110
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <- 0:110
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

compute_fitted_crossedsensis_Infant<- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period2)])), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period1)])), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period1)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27, as.character(period1)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <- 0:110
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <- 0:110
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}


compute_delta_Ex_crossedsensis <- function (period1 = 2000, period2= 2017){
  
  Fitted_sensis_listHump <- compute_fitted_crossedsensis_Hump(period1, period2)
  Fitted_sensis_listSenescent <- compute_fitted_crossedsensis_Senescent(period1, period2)
  Fitted_sensis_listInfant <- compute_fitted_crossedsensis_Infant(period1, period2)
  
  Fitted_DataFrame_base <- Fitted_sensis_listHump[['df_base']]
  Fitted_DataFrame_SensisHump <- Fitted_sensis_listHump[['df_sensis']]
  Fitted_DataFrame_SensisSenescent <- Fitted_sensis_listSenescent[['df_sensis']]
  Fitted_DataFrame_SensisInfant <- Fitted_sensis_listInfant[['df_sensis']]
  
  
  LE_base <- LE_period_Model(QxModel = SSE_deathrates_female_df,Age =  0, Period = period1)
  LE_after <- LE_period_Model(QxModel = SSE_deathrates_female_df,Age =  0, Period = period2)
  
  LE_sensisHump <- LE_period_Model(QxModel = Fitted_DataFrame_SensisHump,Age =  0, Period = 'FittedCurve')
  LE_sensisSenescent <- LE_period_Model(QxModel = Fitted_DataFrame_SensisSenescent,Age =  0, Period = 'FittedCurve')
  LE_sensisInfant <- LE_period_Model(QxModel = Fitted_DataFrame_SensisInfant,Age =  0, Period = 'FittedCurve')
  
  LE_delta_crossed <- as.data.frame(rbind(LE_sensisHump, LE_sensisSenescent, LE_sensisInfant))
  
  
  cbind( LE_base, t(LE_delta_crossed), LE_after)
}


plot_qx_crossedsensis <- function (period1 = 2000, period2 = 2017, component){
  if (component == 'hump'){
    Fitted_sensis_list <- compute_fitted_crossedsensis_Hump(period1, period2)
  }
  else if (component == 'infant'){
    Fitted_sensis_list <- compute_fitted_crossedsensis_Infant(period1, period2)
  }
  else{
    Fitted_sensis_list <- compute_fitted_crossedsensis_Senescent(period1, period2)
  }
  
  
  Fitted_DataFrame_Sensis1 <- Fitted_sensis_list[['df_sensis']]
  Fitted_DataFrame <- Fitted_sensis_list[['df_base']]
  
  Fitted_DataFrame_plot_sensis1 <- melt(Fitted_DataFrame_Sensis1, id.vars = c('Sensis', 'Age'), variable.name = 'Component', value.name = 'qx')
  Fitted_DataFrame_plot <- melt(Fitted_DataFrame, id.vars = c('Sensis','Age'), variable.name = 'Component', value.name = 'qx')
  
  Final_DF_plot <- rbind(Fitted_DataFrame_plot, Fitted_DataFrame_plot_sensis1)
  
  p<- ggplot(Final_DF_plot) + geom_line(aes(x = Age, y = qx, group = interaction(Sensis, Component), color = Component, linetype = Sensis)) +
    scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ') +
    ggtitle(paste0('Components, Data and Fitted Curve for France ', period1, ' - SSE'))
  
  p
}