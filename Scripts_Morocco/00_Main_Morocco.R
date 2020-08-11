######## SCRIPT MODEL SSE data for France #######-
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


####################################################################################-
#---------------------------2. FIT  ----------------------------
####################################################################################-

source('Scripts_Morocco/00_DataPrep_Morocco.R')
source('Scripts_Morocco/01_Model_Fit_Morocco.R')

####################################################################################-
#---------------------------3. SPLINES / COEFF PLOTS  ----------------------------
####################################################################################-
country_code <-  'MAR'
source('Scripts_Morocco/02_Plots_Morocco.R')


####################################################################################-
#---------------------------5. COMPONENTS COMPUTATION / PLOTS  ----------------------------
####################################################################################-
source('Scripts_Morocco/05_Functions_Plots_Morocco.R')
source('Scripts_Morocco/05_Functions_Plots_Females_Morocco.R')

#Plot all years of each component
plot_allyears_senescent()+  scale_color_gradientn(colours = rainbow(5))
plot_allyears_hump() +  scale_color_gradientn(colours = rainbow(5))
plot_allyears_infant()+  scale_color_gradientn(colours = rainbow(5))


#H / M qx plots
F_qx_plot <- plot_allyears_qx() + ggtitle('Female rates')
M_qx_plot <- plot_allyears_qx() + ggtitle('Male rates')
png('MAR_fit_allmx.png', width= 3000, height=2500, res=300)
grid.arrange(F_qx_plot, M_qx_plot)
dev.off()




# Plot all years all components all sex on ONE PAGE
source('Scripts_Morocco/05_Functions_Plots_Females_Morocco.R')
female_infant <- plot_allyears_infant() + ggtitle('Female') + theme(legend.position = 'none')
female_senescent <- plot_allyears_senescent() + theme(legend.position = 'none')
female_hump <- plot_allyears_hump()#+ theme(legend.position = 'none')
gridFemale <- grid.arrange(female_infant, female_senescent, female_hump, ncol=1)

source('Scripts_Morocco/05_Functions_Plots_Morocco.R')
male_infant <- plot_allyears_infant() + ggtitle('Male') + theme(legend.position = 'none')
male_senescent <- plot_allyears_senescent() 
male_hump <- plot_allyears_hump()+ theme(legend.position = 'none')
gridMale <- grid.arrange(male_infant, male_senescent, male_hump, ncol=1)


png('allcompallyearsMar.png', width= 2600, height=3000, res=250)
grid.arrange(gridFemale, gridMale, nrow =1)
dev.off()




png('plotcomponentsyears2.png', width= 3500, height=2000, res=400)
#test PLOT SEVERAL YAERS ALL COMPONENTS
plot_selectedyears_allComp(c(2000,2001))
dev.off()


## All plots
grid.arrange(plot_allyears_senescent(), plot_allyears_hump(),plot_allyears_infant(), plot_allyears_qx())

#SUrvival function plot
plot_allyears_lx()
plot_allyears_qx()


source('Scripts_Morocco/05_Functions_Plots_Females_Morocco.R')
LX_F <- plot_selectyears_lx(c(2000,2004,2016)) + ggtitle('Femmes') + theme(legend.position = 'none')
 
source('Scripts_Morocco/05_Functions_Plots_Morocco.R')
LX_M <- plot_selectyears_lx(c(2000,2004,2016)) + ggtitle('Hommes')

png(file = 'lx_plot.png', width = 13, height = 4, units = 'in', res = 600)
grid.arrange(LX_F, LX_M, nrow=1) 
dev.off()


####################################################################################-
#---------------------------8. COMPONENTS SENSITIVITIES  ----------------------------
####################################################################################-
#Run of sensitivities for EACH components 
 -
#### Plot of components AND raw data AND Fitted Curve
plot_allcompRawData(2016)

#ggsave('plot_allfittedComp.pdf', plot_allFittedcomp, height = 6, width =7 )


#### SENSITIVITY ON COMPONENTS

source('Scripts_Morocco/04_SensisFunctions_Morocco.R')
source('Scripts_Morocco/04_SensisFunctions_Females_Morocco.R')


#Sensitivités -

#Fonctions tests
compute_delta_entropy_sensis(sensis_infant = 0.1, period = 2000)
compute_delta_entropy_sensis(sensis_infant = 10, period = 2016)
compute_delta_Ex_sensis(sensis_infant = 100,period = 2016)
plot_lx_sensis(sensis_infant =  100, sensis_hump = 0.1, sensis_senescent = 0)
plot_qx_sensis(sensis_infant =  100, period = 2015)



#Espérance de vie en stock avec modèle de lissage

source('Scripts_Morocco/04_SensisFunctions_Morocco.R')
sensis_ex_plotdf_m <- plot_Ex_stock()
source('Scripts_Morocco/04_SensisFunctions_Females_Morocco.R')
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



# Crossed sensitivities --

#Plot crossed sensis
source('Scripts_Morocco/04_SensisFunctions_Morocco.R')
sensis_crossed_plot_m  <- plot_Ex_by_comp()
  source('Scripts_Morocco/04_SensisFunctions_Females_Morocco.R')
sensis_crossed_plot_f  <- plot_Ex_by_comp()

sensis_crossed_plot <- rbind(sensis_crossed_plot_f, sensis_crossed_plot_m)

p2 <- ggplot(sensis_crossed_plot) + geom_line(aes(x= year, y= value, group = variable, color= variable, linetype= variable)) +  facet_wrap(~ sex, nrow=2)+
  ggtitle('Component Yearly Ex Improvements')+ ylab('LE improvement (months)') +
  scale_linetype_manual('LE Improvement',values = c(1,1,1,2), labels = c('Hump', 'Senescent', 'Infant', 'Total'))+
  scale_color_discrete('LE Improvement',  labels = c('Hump', 'Senescent', 'Infant', 'Total'))
p2

png('comp_imp_MAR.png', width = 5000, height = 4000, res = 500)
p2
dev.off()




#write.csv(sensis_crossed_df, 'sensis_ex_crossed.csv')
#ggplot(sensis_crossed_df) + geom_line(aes(x= year, y= val)) #+ facet_wrap(~ variable)

png('comp_evolution_graph_MAR.png',  width = 5000, height = 3000, res = 600)
plot_selectedyears_allComp(c(2000,2016)) + ggtitle('Composantes SSE Maroc en 2000 et 2016') + 
  scale_color_discrete('Component', labels=c('Infant', 'Senescent', 'Hump', 'All'))
dev.off()

plot_qx_crossedsensis(2000,2016,component = 'senescent')


imp_corr <- rbind(melt(cor(dcast(subset(sensis_crossed_plot), formula = year ~ variable, value.var = 'value')[,c(2,3,4)])),
            melt(cor(dcast(subset(sensis_crossed_plot), formula = year ~ variable, value.var = 'value')[,c(2,3,4)])))
imp_corr[1:9, 'sex'] <- 'F'
imp_corr[10:18, 'sex'] <- 'M'


corr_F <- ggplot(data = imp_corr, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()  +
  scale_fill_gradient2(midpoint = 0, low = "#52BE80", mid = "white",
                       high = "#03A9F4", space = "Lab" )+
  scale_y_discrete('',labels=c('Hump improvements','Senescent improvements','Infant improvements'))+
  scale_x_discrete('',labels=c('Hump improvements','Senescent improvements','Infant improvements'))+
  geom_text(aes(Var2, Var1, label = round(value,3)), color = "#5D6D7E", size = 5) +
  facet_wrap(~ sex)


png('corrplot_MAR.png', height=1600, width = 3500, res=250)
corr_F
dev.off()



####################################################################################-
#---------------------------9. MODEL ROBUSTESS   ----------------------------
####################################################################################-



