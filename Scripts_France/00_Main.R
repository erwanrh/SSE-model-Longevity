######## SCRIPT MODEL SSE data for France #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# CaFRAda (2016)


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


####################################################################################-
#---------------------------2. FIT  ----------------------------
####################################################################################-

source('Scripts_France/01_Model_Fit.R')

####################################################################################-
#---------------------------3. PLOTS  ----------------------------
####################################################################################-
country_code <-  'FR'
source('Scripts_France/02_Plots.R')


####################################################################################-
#---------------------------4. MUHAT COMPUTATION  ----------------------------
####################################################################################-
source('Scripts_France/05_Functions_Plots.R')
source('Scripts_France/06_Functions_Plots_Females.R')

plot_allyears_senescent()
plot_allyears_hump()
plot_allyears_infant()+  scale_color_gradientn(colours = rainbow(5))

female_infant <- plot_allyears_infant() + ggtitle('Female') + theme(legend.position = 'none')
female_senescent <- plot_allyears_senescent() + theme(legend.position = 'none')
female_hump <- plot_allyears_hump()+ theme(legend.position = 'none')
gridFemale <- grid.arrange(female_infant, female_senescent, female_hump, ncol=1)


male_infant <- plot_allyears_infant() + ggtitle('Male') + theme(legend.position = 'none')
male_senescent <- plot_allyears_senescent() 
male_hump <- plot_allyears_hump()+ theme(legend.position = 'none')
gridMale <- grid.arrange(male_infant, male_senescent, male_hump, ncol=1)


grid.arrange(gridFemale, gridMale, ncol=2)


plot_selectedyears_allComp(c(1965,1966))

grid.arrange(female_qx +  scale_color_gradientn(colours = rainbow(5)),
             male_qx +  scale_color_gradientn(colours = rainbow(5)))
## All plots
grid.arrange(plot_allyears_senescent(), plot_allyears_hump(),plot_allyears_infant(), plot_allyears_qx())

#SUrvival function plot
plot_allyears_lx()
plot_allyears_qx()


source('06_Functions_Plots_Females.R')
LX_F <- plot_selectyears_lx(c(1960,2017)) + ggtitle('Femmes') + theme(legend.position = 'none')+
  annotate("segment", x = 77, xend = 87.5 ,y = 50000, yend = 50000,  size=0.5, arrow=arrow())
source('05_Functions_Plots.R')
LX_M <- plot_selectyears_lx(c(1960,2017)) + ggtitle('Hommes') +
  annotate("segment", x = 70, xend = 82, y = 50000, yend = 50000,  size=0.5, arrow=arrow())

png(file = 'lx_plot.png', width = 13, height = 4, units = 'in', res = 600)
grid.arrange(LX_F, LX_M, nrow=1) 
dev.off()


png('comp_evolution_graphFR.png',  width = 5000, height = 3000, res = 600)
plot_selectedyears_allComp(c(2000,2017)) + ggtitle('Composantes SSE France en 2000 et 2017') + 
  scale_color_discrete('Component', labels=c('Infant', 'Senescent', 'Hump', 'All'))
dev.off()


####################################################################################-
#--------------------------7. SENSITIVITES  ----------------------------
####################################################################################-
#Run of 2 by 2 sensitivities for ALL splines and all coefficients
source('03_Sensis.R')

####################################################################################-
#---------------------------8. COMPONENT SENTIVITIES  ----------------------------
####################################################################################-
#Run of sensitivities for EACH components 
#Run of sensitivities for EACH components 
-
  #### SENSITIVITY ON COMPONENTS
  source('Scripts_France/04_SensisFunctions.R')
source('Scripts_France/04_SensisFunctions_Females.R')

#### Plot of components AND raw data AND Fitted Curve
plot_allcompRawData(2016)


#Sensitivités -

#Fonctions tests
compute_delta_entropy_sensis(sensis_infant = 0.1, period = 2000)
compute_delta_entropy_sensis(sensis_infant = 10, period = 2016)
compute_delta_Ex_sensis(sensis_infant = 100,period = 2016)
plot_qx_sensis(sensis_infant =  100, period = 2015)



#Espérance de vie en stock avec modèle de lissage

source('Scripts/04_SensisFunctions.R')
sensis_ex_plotdf_m <- plot_Ex_stock()
source('Scripts_France/04_SensisFunctions_Females_France.R')
sensis_ex_plotdf_f <- plot_Ex_stock()

sensis_ex_plotdf <- rbind(sensis_ex_plotdf_f,sensis_ex_plotdf_m)

ex_stock_plot <- ggplot(sensis_ex_plotdf,aes(x = year, y= Ex, color= variable),linetype = 1) +
  facet_wrap(~ sex, ncol = 1) +
  geom_point() +
  geom_smooth( method = 'lm', formula = y ~ x) +
  scale_x_continuous(breaks = c(seq(1960, 2015,5),2017)) + 
  ggtitle('Espérance de vie en stock par composante') +
  scale_color_manual(name='Composante',values =  c("#C39BD3", "#5FD69C"), )+
  ylab('Stock LE (months)')

ex_stock_plot 

png('ex_stock.png', height =3000, width = 4000, res = 400)
ex_stock_plot
dev.off()



# Yearly Improvements for each component --

#Plot crossed sensis
source('Scripts_France/04_SensisFunctions.R')
sensis_crossed_plot_m  <- plot_Ex_by_comp()
source('Scripts_France/04_SensisFunctions_Females.R')
sensis_crossed_plot_f  <- plot_Ex_by_comp()

sensis_crossed_plot <- rbind(sensis_crossed_plot_f, sensis_crossed_plot_m)

p2 <- ggplot(subset(sensis_crossed_plot)) + geom_line(aes(x= year, y= value, group = variable, color= variable, linetype= variable)) +  facet_wrap(~ sex, nrow=2)+
  ggtitle('Component Yearly Ex Improvements')+ ylab('LE improvement (months)') +
  scale_linetype_manual('LE Improvement',values = c(1,1,1,2), labels = c('Hump', 'Senescent', 'Infant', 'Total'))+
  scale_color_discrete('LE Improvement',  labels = c('Hump', 'Senescent', 'Infant', 'Total'))
p2

png('comp_imp_FRA4.png', width = 4000, height = 4000, res = 500)
p2
dev.off()


#Evolution of components and fits 
png('comp_evolution_graph_FRA.png',  width = 5000, height = 3000, res = 600)
plot_selectedyears_allComp(c(2000,2016)) + ggtitle('Composantes SSE France en 2000 et 2016') + 
  scale_color_discrete('Component', labels=c('Infant', 'Senescent', 'Hump', 'All'))
dev.off()

plot_qx_crossedsensis(2000,2016,component = 'senescent')



#Corrélatins entre améliorations
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

ggsave(filename = 'corrplot_FRA1.png', plot = corr_F, height=10, width = 7, bg = 'transparent')




####################################################################################-
#---------------------------9. MODEL ROBUSTESS   ----------------------------
####################################################################################-
source('Scripts_France//07_Robustess_France.R')


####################################################################################-
#---------------------------10. OPTIMIZATION   ----------------------------
####################################################################################-
source('Scripts_France//08_OptimFun_France.R')

p1 <- plot_allcompRawData(2016)
ggsave(p1, filename='2016_allcomp_plot.png', width = 8, height = 5)

optimization1 <- optimize_LE(target_LE = 7, perc_Infant = 0.9)
compute_delta_LE_months(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , period= 2016)


p2 <- plot_qx_sensis(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , period= 2016)
ggsave(p2, filename='2016_optimization.png', width = 8, height = 5)

qx_base <- compute_fitted_sensis(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , period= 2016)$df_base$FittedCurve
qx_proj <- compute_fitted_sensis(sensis_infant = optimization1$Infant_imp, sensis_hump = optimization1$Hump_imp , period= 2016)$df_sensis$FittedCurve

RR <- qx_proj / qx_base

