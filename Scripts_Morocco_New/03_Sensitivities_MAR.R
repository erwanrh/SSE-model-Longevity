######## SCRIPT MODEL SSE Sensitivities Analysis #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)




####################################################################################-
#-------------------------- FUNCTIONS SENSITIVITES  ----------------------------
####################################################################################-

# Crossed Sensitivities Males ---------------------------------------------------

compute_fitted_crossedsensis_Senescent_M <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period1)])), 
                                               exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period2)])), 
                                               exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period1)])), 
                                       exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period1)])), 
                                       exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <-0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <-0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

compute_fitted_crossedsensis_Hump_M <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period1)])), 
                                               exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period1)])), 
                                               exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period2)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period1)])), 
                                       exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period1)])), 
                                       exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <-0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <-0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

compute_fitted_crossedsensis_Infant_M <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period2)])), 
                                               exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period1)])), 
                                               exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period1)])), 
                                       exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period1)])), 
                                       exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <-0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <-0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

compute_delta_Ex_crossedsensis_M <- function (period1 = 2000, period2= 2017){
  
  Fitted_sensis_listHump <- compute_fitted_crossedsensis_Hump_M(period1, period2)
  Fitted_sensis_listSenescent <- compute_fitted_crossedsensis_Senescent_M(period1, period2)
  Fitted_sensis_listInfant <- compute_fitted_crossedsensis_Infant_M(period1, period2)
  
  Fitted_DataFrame_base <- Fitted_sensis_listHump[['df_base']]
  Fitted_DataFrame_SensisHump <- Fitted_sensis_listHump[['df_sensis']]
  Fitted_DataFrame_SensisSenescent <- Fitted_sensis_listSenescent[['df_sensis']]
  Fitted_DataFrame_SensisInfant <- Fitted_sensis_listInfant[['df_sensis']]
  
  
  LE_base <- LE_period_Model(QxModel = SSE_deathrates_male_df,Age =  0, Period = period1)
  LE_after <- LE_period_Model(QxModel = SSE_deathrates_male_df,Age =  0, Period = period2)
  
  LE_sensisHump <- LE_period_Model(QxModel = Fitted_DataFrame_SensisHump,Age =  0, Period = 'FittedCurve')
  LE_sensisSenescent <- LE_period_Model(QxModel = Fitted_DataFrame_SensisSenescent,Age =  0, Period = 'FittedCurve')
  LE_sensisInfant <- LE_period_Model(QxModel = Fitted_DataFrame_SensisInfant,Age =  0, Period = 'FittedCurve')
  
  LE_delta_crossed <- as.data.frame(rbind(LE_sensisHump, LE_sensisSenescent, LE_sensisInfant))
  
  
  cbind( LE_base, t(LE_delta_crossed), LE_after)
}

# Crossed Sensitivities Females ---------------------------------------------------

compute_fitted_crossedsensis_Senescent_F <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period2)])), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period1)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <-0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <-0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

compute_fitted_crossedsensis_Hump_F <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period1)])), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period2)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period1)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <-0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <-0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

compute_fitted_crossedsensis_Infant_F <- function( period1 = 2000, period2 = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period2)])), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period1)])), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period1)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period1)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period1)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <-0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <-0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}


compute_delta_Ex_crossedsensis_F <- function (period1 = 2000, period2= 2017){
  
  Fitted_sensis_listHump <- compute_fitted_crossedsensis_Hump_F(period1, period2)
  Fitted_sensis_listSenescent <- compute_fitted_crossedsensis_Senescent_F(period1, period2)
  Fitted_sensis_listInfant <- compute_fitted_crossedsensis_Infant_F(period1, period2)
  
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

# Sensitivities Males -----------

compute_fitted_sensis_M <- function(sensis_infant=0, sensis_hump=0, sensis_senescent=0, period = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period)])*(1+sensis_infant)), 
                                               exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period)])*(1+sensis_senescent)), 
                                               exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period)])*(1+sensis_hump))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2, as.character(period)])), 
                                       exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:22, as.character(period)])), 
                                       exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[23:42, as.character(period)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <- 0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <- 0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}

# Sensitivities Females -----------

compute_fitted_sensis_F <- function(sensis_infant=0, sensis_hump=0, sensis_senescent=0, period = 2017){
  
  Fitted_DataFrame_Sensis1 <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period)])*(1+sensis_infant)), 
                                               exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period)])*(1+sensis_senescent)), 
                                               exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period)])*(1+sensis_hump))))
  colnames(Fitted_DataFrame_Sensis1) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame <- data.frame(cbind(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2, as.character(period)])), 
                                       exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:22, as.character(period)])), 
                                       exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[23:42, as.character(period)]))))
  colnames(Fitted_DataFrame) <- c('Infant', 'Senescent', 'Hump')
  
  Fitted_DataFrame_Sensis1$FittedCurve <- rowSums(Fitted_DataFrame_Sensis1)
  Fitted_DataFrame_Sensis1$Age <- 0:85
  Fitted_DataFrame_Sensis1$Sensis <- 'Sensis'
  
  Fitted_DataFrame$FittedCurve <- rowSums(Fitted_DataFrame)
  Fitted_DataFrame$Age <- 0:85
  Fitted_DataFrame$Sensis <- 'Base'
  
  list('df_base'= Fitted_DataFrame, 'df_sensis'= Fitted_DataFrame_Sensis1)
}


# PLOT ---------------------
compute_Ex_by_comp <- function() {
  
  
  sensis_crossed_df_M<- data.frame()
  sensis_crossed_df_F <- data.frame()
  
  year_1 <- as.numeric(colnames(SSE_deathrates_male_df)[1])
  
  for (year in as.numeric(colnames(SSE_deathrates_male_df)[-1]))
  {
    sensis_crossed_df_M <- rbind(sensis_crossed_df_M,cbind(compute_delta_Ex_crossedsensis_M(period1 = year_1, period2 = year), year))
    sensis_crossed_df_F <- rbind(sensis_crossed_df_F,cbind(compute_delta_Ex_crossedsensis_F(period1 = year_1, period2 = year), year))
    year_1 <- year
  }
  
  
  sensis_crossed_var_df_M <- (sensis_crossed_df_M[, c('LE_sensisHump','LE_sensisSenescent','LE_sensisInfant')] - sensis_crossed_df_M$LE_base)*12
  sensis_crossed_var_df_M$year <- sensis_crossed_df_M$year
  sensis_crossed_var_df_M$LE_var <- (sensis_crossed_df_M$LE_after - sensis_crossed_df_M$LE_base)*12
  
  sensis_crossed_var_df_F <- (sensis_crossed_df_F[, c('LE_sensisHump','LE_sensisSenescent','LE_sensisInfant')] - sensis_crossed_df_F$LE_base)*12
  sensis_crossed_var_df_F$year <- sensis_crossed_df_F$year
  sensis_crossed_var_df_F$LE_var <- (sensis_crossed_df_F$LE_after - sensis_crossed_df_F$LE_base)*12
  
  
  sensis_crossed_plot_m  <- melt(sensis_crossed_var_df_M, id.vars = 'year' )
  sensis_crossed_plot_m$sex <- 'M'
  
  sensis_crossed_plot_f  <- melt(sensis_crossed_var_df_F, id.vars = 'year' )
  sensis_crossed_plot_f$sex <- 'F'
  
  rbind(sensis_crossed_plot_m, sensis_crossed_plot_f)
}


plot_Ex_by_comp<- function(){
  
  sensis_crossed_plot <- compute_Ex_by_comp()
  
  p2 <- ggplot(subset(sensis_crossed_plot)) + geom_line(aes(x= year, y= value, group = variable, color= variable, linetype= variable)) +  facet_wrap(~ sex, nrow=2)+
    ggtitle('Component Yearly Ex Improvements')+ ylab('LE improvement (months)') +
    scale_linetype_manual('LE Improvement',values = c(1,1,1,2), labels = c('Hump', 'Senescent', 'Infant', 'Total'))+
    scale_color_discrete('LE Improvement',  labels = c('Hump', 'Senescent', 'Infant', 'Total'))
  
  p2
  
}

compute_Ex_stock <- function(){
  #Dataframe des sensibilités d'ex
  sensis_ex_dfF <- data.frame()
  sensis_ex_dfM <- data.frame()
  
  for(delta in c(100)){
    for (year in 2000:2016){
      sensis_ex_dfF[as.character(paste(year,delta, sep = '_')), 'Infant'] <- compute_delta_Ex_sensis_F(sensis_infant = delta, period = year)*12
      sensis_ex_dfF[as.character(paste(year,delta, sep = '_')), 'Hump'] <- compute_delta_Ex_sensis_F(sensis_hump = delta, period = year)*12
      sensis_ex_dfF[as.character(paste(year,delta, sep = '_')), 'Senescent'] <- compute_delta_Ex_sensis_F(sensis_senescent = delta, period = year)*12
      sensis_ex_dfF[as.character(paste(year,delta, sep = '_')), 'variation'] <- delta
      sensis_ex_dfF[as.character(paste(year,delta, sep = '_')), 'year'] <- year
      
      sensis_ex_dfM[as.character(paste(year,delta, sep = '_')), 'Infant'] <- compute_delta_Ex_sensis_M(sensis_infant = delta, period = year)*12
      sensis_ex_dfM[as.character(paste(year,delta, sep = '_')), 'Hump'] <- compute_delta_Ex_sensis_M(sensis_hump = delta, period = year)*12
      sensis_ex_dfM[as.character(paste(year,delta, sep = '_')), 'Senescent'] <- compute_delta_Ex_sensis_M(sensis_senescent = delta, period = year)*12
      sensis_ex_dfM[as.character(paste(year,delta, sep = '_')), 'variation'] <- delta
      sensis_ex_dfM[as.character(paste(year,delta, sep = '_')), 'year'] <- year
    }
    
  }
  
  #Plot des sensibilités d'ex pour les hommes et les femmes 
  sensis_ex_plotdf_m <- melt(sensis_ex_dfM, id.vars = 'year', measure.vars = c('Infant', 'Hump'), value.name = 'Ex')
  sensis_ex_plotdf_m$sex = 'M'
  
  sensis_ex_plotdf_f <- melt(sensis_ex_dfF, id.vars = 'year', measure.vars = c('Infant', 'Hump'), value.name = 'Ex')
  sensis_ex_plotdf_f$sex = 'F'
  
  rbind(sensis_ex_plotdf_f, sensis_ex_plotdf_m) 
}



plot_Ex_stock <- function(){
  sensis_ex_plotdf <- compute_Ex_stock()
  
  ex_stock_plot <- ggplot(sensis_ex_plotdf, aes(x = year, y = Ex, color = variable), linetype = 1) +
    facet_wrap(~ sex, ncol = 1) +
    geom_point() +
    geom_smooth( method = 'lm', formula = y ~ x) +
    scale_x_continuous(breaks = c(seq(1960, 2015,5),2017)) + 
    ggtitle('Espérance de vie en stock par composante') +
    scale_color_manual(name='Composante',values =  c("#C39BD3", "#5FD69C"), )+
    ylab('Stock LE (months)')
  
  ex_stock_plot 
}

# CORRELATIONS ---------------------

plot_ExImp_corr <- function(){
  #Corrélatins entre améliorations
  sensis_crossed_plot <- compute_Ex_by_comp()
  
  imp_corr <- rbind(melt(cor(dcast(subset(sensis_crossed_plot, sex == 'F'), formula = year ~ variable, value.var = 'value')[,c(2,3,4)])),
                    melt(cor(dcast(subset(sensis_crossed_plot, sex == 'M'), formula = year ~ variable, value.var = 'value')[,c(2,3,4)])))
  imp_corr[1:9, 'sex'] <- 'F'
  imp_corr[10:18, 'sex'] <- 'M'
  
  
  corr_F <- ggplot(data = imp_corr, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()  +
    scale_fill_gradient2(midpoint = 0, low = "#52BE80", mid = "white",
                         high = "#03A9F4", space = "Lab" )+
    scale_y_discrete('',labels=c('Hump','Senescent','Infant'))+
    scale_x_discrete('',labels=c('Hump','Senescent','Infant'))+
    geom_text(aes(Var2, Var1, label = round(value,3)), color = "#5D6D7E", size = 5) +
    facet_wrap(~ sex, nrow = 2) + theme(rect = element_rect(fill = "transparent") # all rectangles
    )
  corr_F
  
}




