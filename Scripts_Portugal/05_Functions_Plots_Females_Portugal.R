######## SCRIPT MODEL SSE   #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)




####################################################################################-
#-------------------------- FUNCTIONS TO PLOT FOR FEMALES  ----------------------------
####################################################################################-

plot_allyears_senescent <- function(){
  #PLot du comportement de la mortalité sénescente
  senescent <- melt(exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescent$Comp <- '2'
  ggplot(senescent) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 
}

plot_allyears_hump <- function(){
  #PLot du comportement de la bosse
  hump <- melt(exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  hump$Comp <- '3'
  ggplot(hump) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10', limits = c(1e-06, 1)) 
}

plot_allyears_infant <- function(){
  #PLot du comportement de la mortalité infantile
  infant_qx <- exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2,]))
  
  first_infant_rate <- infant_qx[1,]
  
  infant <- melt(infant_qx,varnames = c('Age','Year'),value.name = 'm')
  infant$Comp <- '1'
  ggplot(infant) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') #+
   # annotation_custom(tableGrob(round(first_infant_rate, 4),rows = names(first_infant_rate), cols = ('q0') , theme=ttheme_default(base_size = 7)), 
    #                  xmin=80, xmax=110)
  
}

plot_allyears_qx <- function(){
  hump <- melt(exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  hump$Comp <- '3'
  senescent <- melt(exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescent$Comp <- '2'
  #PLot du comportement de la mortalité infantile
  infant_qx <- exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2,]))
  
  first_infant_rate <- infant_qx[1,]
  
  infant <- melt(infant_qx,varnames = c('Age','Year'),value.name = 'm')
  infant$Comp <- '1'
  #PLot du comportement des Qx Fittés
  MU <- hump
  MU$m <- hump$m+infant$m+senescent$m
  MU$Comp <-'All'
  ggplot(MU) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') +
    scale_color_gradientn(colours = rainbow(5))
}

plot_selectedyears_allComp <- function(Years){
  hump <- melt(exp(SSE_females$XX$X3%*%(SSE_coeffcients_females_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
  hump$Comp <- '3'
  senescent <- melt(exp(SSE_females$XX$X2%*%(SSE_coeffcients_females_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
  senescent$Comp <- '2'
  infant_qx <- exp(SSE_females$XX$X1%*%(SSE_coeffcients_females_df[1:2,]))
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
  
  Fitted_DataFrame <- compute_fitted_sensis(period = period)[['df_base']]
  
  Fitted_DataFrame_plot <- melt(Fitted_DataFrame[,-6], id.vars = 'Age', variable.name = 'Component', value.name = 'qx')
  Fitted_DataFrame_plot$Data <- SSE_data_male_df[, as.character(period)]
  Fitted_DataFrame_plot$Data1 <- 'Data'
  
  plot_allFittedcomp <- ggplot(Fitted_DataFrame_plot) + geom_line(aes(x = Age, y = qx, group = Component, color = Component)) +
    geom_point(aes(x= Age, y = Data), shape = 20,  size = 3, color = 'purple', alpha = 0.1)+
    scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ') +
    ggtitle(paste0('Components, Data and Fitted Curve for France ', period,' - SSE'))
  plot_allFittedcomp
  
}

plot_allyears_lx <- function(){
  row.names(SSE_deathrates_female_df) <- 0:85
  
  lx <- compute_lx_period(SSE_deathrates_female_df, start.age = 0)
  lx_plot <- melt(lx, varnames = c('Age', 'Period'), value.name = 'lx')
  
  entropy <- compute_entropy_period(SSE_deathrates_female_df, start.age = 0)
  
  
  p <- ggplot(lx_plot) + geom_line(aes(x=Age, y = lx, group = Period, color = Period)) +
    annotation_custom(tableGrob(cbind(as.character(entropy[,1]), round(entropy[,2]*100,3)), 
                                rows = NULL, theme=ttheme_default(base_size = 7), cols = c('Period', 'Entropy')),
                      xmin=0, xmax=30, ymin=20000, ymax=25000) +
    ggtitle('Survival curves and entropy on fitted values')
  
  p 
}


plot_selectyears_lx <- function(years){
  row.names(SSE_deathrates_female_df) <- 0:85
  
  lx <- compute_lx_period(SSE_deathrates_female_df, start.age = 0)
  lx_plot <- melt(lx, varnames = c('Age', 'Period'), value.name = 'lx')
  
  entropy <- compute_entropy_period(SSE_deathrates_female_df[,colnames(SSE_deathrates_female_df) %in% years], start.age = 0)
  
  
  p <- ggplot(subset(lx_plot, Period %in% years)) + geom_line(aes(x=Age, y = lx, group = as.character(Period), color = as.character(Period))) +
    annotation_custom(tableGrob(cbind(as.character(entropy[,1]), round(entropy[,2]*100,3)), 
                                rows = NULL, theme=ttheme_default(base_size = 7), cols = c('Period', 'Entropy')),
                      xmin=0, xmax=30, ymin=20000, ymax=25000) +
    ggtitle('Survival curves and entropy on fitted values')
  
  p 
}


