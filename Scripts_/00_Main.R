######## SCRIPT SSE #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)


####################################################################################-
#---------------------------1. INITIALIZATION  ----------------------------
####################################################################################-
rm(list=ls()); gc(); #Clean environment

#install.packages('HMDHFDplus')
library(HMDHFDplus)
library(Matrix)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(zoo)         

# Mortality Smoothing Package (required ro run Mort Hump functions 
pckg<-read.csv('MortalitySmooth-master/pckSSE.csv')
for (pk in pckg[,1]){
  source(paste0('MortalitySmooth-master/MortalitySmooth-master/R/',pk,'.R'))  
}

# Mort Hump package 
source('MortalitySmooth-master/SSE_fit.R')
#Indicators functions
source('00_Functions_Calc.R')

#FILL THE RIGHT COUNTRY CODE
country_code <-  'FR'

####################################################################################-
#---------------------------2. FIT  ----------------------------
####################################################################################-
#Fit of model and parameters are located in the following script  
source('Scripts_/01_Model_Fit.R')

####################################################################################-
#---------------------------3. MODEL PLOTS   ----------------------------
####################################################################################-
source('Scripts_/02_Plots.R')

plot_AlphasFitted()
plot_rolling5Y_Ex()

####################################################################################-
#---------------------------4. COMPONENTS PLOTS  ----------------------------
####################################################################################-


plot1 <- plot_allyears_allcompgrid()
ggsave(paste0('allcomponents_', country_code,'.png'), plot = plot1, width = 6, height = 6)

plot_selectedyears_allComp(c(1965,1966))


## All plots
grid.arrange(plot_allyears_senescent(), plot_allyears_hump(),plot_allyears_infant(), plot_allyears_qx())

#SUrvival function plot
plot_allyears_lx()

#QX plots
plot2 <- plot_allyears_qx()
ggsave(paste0('allyearsQx_', country_code,'.png'), plot = plot2, width = 6, height = 6)

#### Plot of components AND raw data AND Fitted Curve
plot_allcompRawData(2010)

####################################################################################-
#--------------------------5. SENSITIVITES  ----------------------------
####################################################################################-
#Run of 2 by 2 sensitivities for ALL splines and all coefficients
source('Scripts_/03_Sensitivities.R')


# Yearly Improvements for each component --
p2 <- plot_Ex_by_comp()

png(paste0('comp_imp_',country_code, 'A5.png'), width = 4000, height = 4000, res = 500)
p2
dev.off()


#Evolution of components and fits 
png(paste0('comp_evolution_graph_',country_code, '.png'),  width = 5000, height = 3000, res = 600)
plot_selectedyears_allComp(c(2000,2016)) + ggtitle('Composantes SSE France en 2000 et 2016') + 
  scale_color_discrete('Component', labels=c('Infant', 'Senescent', 'Hump', 'All'))
dev.off()

plot_qx_crossedsensis(2000,2016,component = 'senescent')


#Corrélation entre améliorations par composantes
corr_plt <- plot_impLEcorr()
ggsave(filename = paste0('imp_correlation_',country_code, '.png'), plot = corr_plt, height=10, width = 7, bg = 'transparent')


#Espérance de vie en stock avec modèle de lissage
plot_Ex_stock()



####################################################################################-
#---------------------------6. OPTIMIZATION   ----------------------------
####################################################################################-
source('Scripts_/04_Optimization.R')

p1 <- plot_allcompRawData(2016)
ggsave(p1, filename=paste0('all_components_rawdata_',country_code, '.png'), width = 8, height = 5)

optimization1 <- optimize_LE(target_LE = 4, perc_Infant = 0.9, SSE_coefficients = SSE_coeffcients_females_df[,'2017'])
compute_delta_Ex_sensis_optim(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp, SSE_coefficients = SSE_coeffcients_females_df[,'2017'])


plot_Ex()
plot_ExImp()


p2 <- plot_qx_sensis(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , period= 2016)
ggsave(p2, filename=paste0('optimization_',country_code, '.png'), width = 8, height = 5)

qx_base <- compute_fitted_sensis(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , period= 2016)$df_base$FittedCurve
qx_proj <- compute_fitted_sensis(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , period= 2016)$df_sensis$FittedCurve

RR <- qx_proj / qx_base

