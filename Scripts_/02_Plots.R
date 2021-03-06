######## SCRIPT MODEL SSE Plot Functions #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)




####################################################################################-
#-------------------------- Plot funtions  ----------------------------
####################################################################################-



plot_AlphasFitted <- function() {
  
  COMP1 <- paste0('X',1:2)
  COMP2 <- paste0('X',3:27)
  COMP3 <- paste0('X',28:52)
  
  Coeff_plot_males <- melt(SSE_coeffcients_males_df,varnames = c('Coeff.num','Year'),value.name = 'coeff',)
  Coeff_plot_females <- melt(SSE_coeffcients_females_df,varnames = c('Coeff.num','Year'),value.name = 'coeff',)
  Coeff_plot_males$gender <- 'M'
  Coeff_plot_females$gender <- 'F'
  Coeff_plot <- rbind(Coeff_plot_males,Coeff_plot_females)
  
  Coeff_plot [Coeff_plot$Coeff.num %in% COMP1, 'comp']<-1
  Coeff_plot[Coeff_plot$Coeff.num %in% COMP2, 'comp']<-2
  Coeff_plot[Coeff_plot$Coeff.num %in% COMP3, 'comp']<-3
  
  coefplots <- ggplot(Coeff_plot)+geom_line(aes(x=Year, y=coeff, color=Coeff.num))+ ggtitle(paste0('SSE Coefficients - ',country_code)) +facet_wrap(gender~comp,scales = 'free')

  coefplots
}

#To save the plot
#ggsave(paste0('coefplots_',country_code,'.pdf'), Plot_AlphasFitted(), width = 20, height = 12)




plot_BSplines <- function(){
  
  splines_plot <- SSE_splines_male
  splines_plot [splines_plot$splinenb %in% 1:2, 'comp']<-'c1'
  splines_plot[splines_plot$splinenb %in% 3:27, 'comp']<-'c2'
  splines_plot[splines_plot$splinenb %in% 28:52, 'comp']<-'c3'
  
  splines_ggplot <- ggplot(splines_plot,aes(x=age, y=splinevalue, color=year, group=year))+ geom_line() + facet_wrap(.~splinenb+ comp) + ggtitle(paste0('Splines - ', country_code))

  splines_ggplot
  
}

#To save plot
#ggsave(paste0('spline_comp_',country_code,'.pdf'),plot = splines_ggplot,width = 15, height = 15)





plot_allyears_senescent <- function(){
  #PLot du comportement de la mortalité sénescente
  senescentM <- melt(exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescentM$Comp <- '2'
  senescentM$sex <- 'M'
  
  senescentF <- melt(exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescentF$Comp <- '2'
  senescentF$sex <- 'F'
  
  senescent <- rbind(senescentM, senescentF)
  
  ggplot(senescent) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') + facet_wrap(sex ~.)
}


plot_allyears_hump <- function(){
  #PLot du comportement de la bosse
  humpM <- melt(exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  humpM$Comp <- '3'
  humpM$sex <- 'M'
  
  humpF <- melt(exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  humpF$Comp <- '3'
  humpF$sex <- 'F'
  
  hump <- rbind(humpM, humpF)
  
  ggplot(hump) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10', limits = c(1e-06, 1)) + facet_wrap( sex ~.) 
}



plot_allyears_infant <- function(){
  #Plot du comportement de la mortalit? infantile
  infantM <- melt(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2,])),varnames = c('Age','Year'),value.name = 'm')
  infantM$Comp <- '1'
  infantM$sex <- 'M'
  
  infantF <- melt(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2,])),varnames = c('Age','Year'),value.name = 'm')
  infantF$Comp <- '1'
  infantF$sex <- 'F'
  
  infant <- rbind(infantM, infantF)
  
  ggplot(infant) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') + facet_wrap(sex ~.)

}



plot_allyears_qx <- function(){
  DeathRates_plot_males <- melt(SSE_deathrates_male_df,varnames = c('Age','Year'),value.name = 'coeff',)
  DeathRates_plot_females <- melt(SSE_deathrates_female_df,varnames = c('Age','Year'),value.name = 'coeff',)
  DeathRates_plot_males$gender <- 'M'
  DeathRates_plot_females$gender <- 'F'
  DeathRates_plot <- rbind(DeathRates_plot_females,DeathRates_plot_males)
  DeathRates_plot$Age <- as.numeric(DeathRates_plot$Age)
  
  DeathRatesplot <- ggplot(DeathRates_plot)+geom_line(aes(x=Age, y=coeff, group = Year,color=Year))+ ggtitle(paste0('SSE Mortality Rates - ', country_code)) +
    facet_wrap(gender~.,scales = 'free', nrow = 2)+
    scale_y_continuous(trans='log2', labels = scales::scientific) +
    scale_x_continuous(breaks = seq(0,110,10), labels= seq(0,110,10) )
  
  DeathRatesplot
}


plot_allyears_allcompgrid <- function(){
  senescentM <- melt(exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescentM$Comp <- 'Senescent'
  senescentM$sex <- 'Male'
  
  senescentF <- melt(exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescentF$Comp <- 'Senescent'
  senescentF$sex <- 'Female'
  
  senescent <- rbind(senescentM, senescentF)
  
  humpM <- melt(exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  humpM$Comp <- 'Hump'
  humpM$sex <- 'Male'
  
  humpF <- melt(exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  humpF$Comp <- 'Hump'
  humpF$sex <- 'Female'
  
  hump <- rbind(humpM, humpF)

  infantM <- melt(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2,])),varnames = c('Age','Year'),value.name = 'm')
  infantM$Comp <- 'Infant'
  infantM$sex <- 'Male'
  
  infantF <- melt(exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2,])),varnames = c('Age','Year'),value.name = 'm')
  infantF$Comp <- 'Infant'
  infantF$sex <- 'Female'
  
  infant <- rbind(infantM, infantF)
  
  all_comp_df <- rbind(senescent, hump, infant)
  
  allcompplot <- ggplot(all_comp_df)+
    geom_line(aes(x=Age, y=m, group = Year, color=Year))+ 
    scale_y_continuous(trans='log2', limits = c(1e-06, 1), labels = scales::scientific) +
    facet_grid(Comp ~sex) + ggtitle(paste0('SSE components for ', country_code))
  
  allcompplot

}




plot_selectedyears_allComp <- function(Years){
  hump <- melt(exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  hump$Comp <- '3'
  senescent <- melt(exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescent$Comp <- '2'
  infant_qx <- exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2,]))
  infant <- melt(infant_qx,varnames = c('Age','Year'),value.name = 'm')
  infant$Comp <- '1'
  #PLot du comportement des Qx Fittés
  MU <- hump
  MU$m <- hump$m+infant$m+senescent$m
  MU$Comp <-'All'
  
  #Plot du comportement de toutes les composantes + des qx fittés
  ALL_plot_df <- rbind(hump, infant, senescent, MU)
  ALL_plot_df$Year <-as.character(ALL_plot_df$Year)
  ggplot(subset(ALL_plot_df, Year %in% Years)) + 
    geom_line( aes(x=Age, y=m, group= interaction(Year,Comp), color=Comp,linetype=Year))+
    scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ')  +ggtitle('SSE Composantes')
}

plot_allcompRawData <- function(period=2017){
  
  Fitted_DataFrameM <- compute_fitted_sensis_M(period = period)[['df_base']]
  
  Fitted_DataFrame_plotM <- melt(Fitted_DataFrameM[,-6], id.vars = 'Age', variable.name = 'Component', value.name = 'qx')
  Fitted_DataFrame_plotM$Data <- SSE_data_male_df[, as.character(period)]
  Fitted_DataFrame_plotM$Data1 <- 'Data'
  Fitted_DataFrame_plotM$sex <- 'Male'
  
  Fitted_DataFrameF <- compute_fitted_sensis_F(period = period)[['df_base']]
  
  Fitted_DataFrame_plotF <- melt(Fitted_DataFrameF[,-6], id.vars = 'Age', variable.name = 'Component', value.name = 'qx')
  Fitted_DataFrame_plotF$Data <- SSE_data_female_df[, as.character(period)]
  Fitted_DataFrame_plotF$Data1 <- 'Data'
  Fitted_DataFrame_plotF$sex <- 'Female'
  
  Fitted_DataFrame_plot <- rbind(Fitted_DataFrame_plotM, Fitted_DataFrame_plotF)
  

  plot_allFittedcomp <- ggplot(Fitted_DataFrame_plot) + geom_line(aes(x = Age, y = qx, group = Component, color = Component)) +
    geom_point(aes(x= Age, y = Data), shape = 20,  size = 3, color = 'purple', alpha = 0.1)+
    scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ') +
    ggtitle(paste0('Components, Data and Fitted Curve - ', country_code, ' / ' ,period,' - SSE')) +
    facet_wrap(. ~ sex, nrow =2 )
  plot_allFittedcomp
  
}

plot_allyears_lx <- function(){
  row.names(SSE_deathrates_male_df) <- 0:110
  
  lx <- compute_lx_period(SSE_deathrates_male_df, start.age = 0)
  lx_plot <- melt(lx, varnames = c('Age', 'Period'), value.name = 'lx')
  
  entropy <- compute_entropy_period(SSE_deathrates_male_df, start.age = 0)
  
  
  p <- ggplot(lx_plot) + geom_line(aes(x=Age, y = lx, group = Period, color = Period)) +
    annotation_custom(tableGrob(cbind(as.character(entropy[,1]), round(entropy[,2]*100,3)), 
                                rows = NULL, theme=ttheme_default(base_size = 7), cols = c('Period', 'Entropy')),
                      xmin=0, xmax=30, ymin=20000, ymax=25000) +
    ggtitle('Survival curves and entropy on fitted values')
  
  p 
}


plot_selectyears_lx <- function(years){
  
  lx <- compute_lx_period(SSE_deathrates_male_df, start.age = 0)
  lx_plot <- melt(lx, varnames = c('Age', 'Period'), value.name = 'lx')
  
  entropy <- compute_entropy_period(SSE_deathrates_male_df[,colnames(SSE_deathrates_male_df) %in% years], start.age = 0)
  
  
  p <- ggplot(subset(lx_plot, Period %in% years)) + geom_line(aes(x=Age, y = lx, group = as.character(Period), color = as.character(Period))) +
    annotation_custom(tableGrob(cbind(as.character(entropy[,1]), round(entropy[,2]*100,3)), 
                                rows = NULL, theme=ttheme_default(base_size = 7), cols = c('Period', 'Entropy')),
                      xmin=0, xmax=30, ymin=20000, ymax=25000) +
    ggtitle('Survival curves and entropy on fitted values') +
    scale_color_discrete(name='Period')
  
  p 
}



# Ex plots -----------

compute_Ex <- function (){
  ex_m <- data.frame()
  ex_f <- data.frame()
  for (year_ in as.numeric(colnames(SSE_deathrates_male_df))){
    ex_m[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_male_df) 
    ex_f[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_female_df)
    
  }
  
  ex_m$var <- c(NA, diff(ex_m$V1))*12
  ex_m$sex <- 'Male'
  ex_m$year <- as.numeric(colnames(SSE_deathrates_male_df))
  
  
  ex_f$var <- c(NA, diff(ex_f$V1))*12
  ex_f$sex <- 'Female'
  ex_f$year <- as.numeric(colnames(SSE_deathrates_male_df))
  
  
  rbind(ex_m, ex_f)
  
}



plot_Ex <- function() {
  ex_tot<- compute_Ex()
  ex_plot <- ggplot(ex_tot) + geom_line(aes(x= year, y=V1, group = sex, color =sex)) + ylab('Ex') +
    theme(legend.position = "none") + ggtitle(paste0('Life expectancy plot - ', country_code))
  ex_plot
  
}


plot_ExImp <- function() {
  ex_tot<- compute_Ex()
  ex_imp_plot <- ggplot(ex_tot, aes(x= year, y = var, group = sex, color = sex)) +
    geom_line() + ylab('LE Improvement (months)') +
    ggtitle('LE Improvements')
  ex_imp_plot
  
}



compute_rolling5Y_Ex <- function(){
  ex_tot <-compute_Ex()
  ex_f <- subset(ex_tot, sex == 'Female')
  ex_m <- subset(ex_tot, sex == 'Male')               
  
  ex_mean_rolling5 <- data.frame()
  i<-1
  for (year in unique(ex_tot$year)){
    ex_mean_rolling5[i,'var'] <- mean(ex_f[ex_f$year <= year & ex_f$year > year-5, 'var'])
    ex_mean_rolling5[i+1,'var'] <- mean(ex_m[ex_m$year <= year & ex_m$year > year-5, 'var'])
    ex_mean_rolling5[i,'year'] <- year
    ex_mean_rolling5[i+1,'year'] <- year
    ex_mean_rolling5[i,'sex'] <- 'Female'
    ex_mean_rolling5[i+1,'sex'] <- 'Male'
    i <- i+2
  }
  ex_mean_rolling5
}


plot_rolling5Y_Ex <- function(){
  ex_mean_rolling5 <- compute_rolling5Y_Ex()
  ex_tot <- compute_Ex()
  ex_imp_plot_5yrolling <- ggplot() +
    geom_line(data = ex_mean_rolling5, aes(x= year, y = var, group = sex,linetype= sex), color = '#5DADE2') +
    geom_point(data= ex_tot, aes(x= year, y = var, group = sex, shape= sex), alpha = 0.5, )+
    scale_shape_manual('Improvements', values = c(16,16)) +
    scale_linetype_manual('Rolling 5Y Average', values = c(1,1)) +
    scale_x_continuous(breaks = seq(min(ex_tot$year),max(ex_tot$year),10), minor_breaks =seq(min(ex_tot$year),max(ex_tot$year),5))+
    ylab('LE Improvement (months)') + 
    ggtitle('LE Improvements and 5Y rolling AVG') +
    facet_wrap(~ sex, nrow = 2)
  
  ex_imp_plot_5yrolling
}




