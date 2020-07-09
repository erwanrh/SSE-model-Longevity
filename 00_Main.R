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

source('01_Model_Fit.R')

####################################################################################-
#---------------------------3. PLOTS  ----------------------------
####################################################################################-
country_code <-  'FR'
source('02_Plots.R')


####################################################################################-
#---------------------------4. MUHAT COMPUTATION  ----------------------------
####################################################################################-
source('05_Functions_Plots.R')
source('06_Functions_Plots_Females.R')
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

grid.arrange(gridFemale, gridMale, nrow =1)
 
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

#### Plot of components AND raw data AND Fitted Curve
plot_allcompRawData(1960)




#### SENSITIVITY ON COMPONENTS

source('04_SensisFunctions.R')
source('04_SensisFunctions_Females.R')


#Sensitivités -------------------

#Fonctions tests
compute_delta_entropy_sensis(sensis_infant = 0.1, period = 1990)
compute_delta_entropy_sensis(sensis_infant = 10, period = 2017)
compute_delta_Ex_sensis(sensis_infant = 100,period = 2017)
plot_lx_sensis(sensis_infant =  100, sensis_hump = 0.1, sensis_senescent = 0)
plot_qx_sensis(sensis_infant =  100, period = 2017)

source('04_SensisFunctions.R')
source('04_SensisFunctions_Females.R')
#Dataframe des sensibilités d'ex
sensis_ex_df <- data.frame()
sensis_entropy_df <- data.frame()

for(delta in c(100)){
  for (year in 1960:2017){
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'Infant'] <- compute_delta_Ex_sensis(sensis_infant = delta, period = year)*12
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'Hump'] <- compute_delta_Ex_sensis(sensis_hump = delta, period = year)*12
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'Senescent'] <- compute_delta_Ex_sensis(sensis_senescent = delta, period = year)*12
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'variation'] <- delta
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'year'] <- year
  }
  
}


write.csv2(sensis_ex_df, 'sensis_life_expectancy2.csv')
write.csv2(sensis_entropy_df, 'sensis_entropy.csv')
#ggsave('plot_allfittedComp_Sensis.pdf', plot_allFittedcomp, height = 6, width =7 )



#Plot des sensibilités d'ex pour les hommes et les femmes 
sensis_ex_plotdf_f <- melt(sensis_ex_df, id.vars = 'year', measure.vars = c('Infant', 'Hump'), value.name = 'Ex')
sensis_ex_plotdf_f$sex = 'F'

sensis_ex_plotdf_m <- melt(sensis_ex_df, id.vars = 'year', measure.vars = c('Infant', 'Hump'), value.name = 'Ex')
sensis_ex_plotdf_m$sex = 'M'

sensis_ex_plotdf <- rbind(sensis_ex_plotdf_f,sensis_ex_plotdf_m)

#Espérance de vie en stock avec modèle de lissage
ex_stock_plot <- ggplot(sensis_ex_plotdf,aes(x = year, y= Ex, color= variable),linetype = 1) +
  facet_wrap(~ sex, ncol = 1) +
  geom_point() +
  geom_smooth( method = 'lm', formula = y ~ x) +
  scale_x_continuous(breaks = c(seq(1960, 2015,5),2017)) + 
  ggtitle('Espérance de vie en stock par composante') +
  scale_color_manual(name='Composante',values =  c("#C39BD3", "#5FD69C"), )+
  ylab('Stock LE (months)')

png('ex_stock.png', height =3000, width = 4000, res = 400)
ex_stock_plot
dev.off()



# Crossed sensitivities ---------------------------------------------------
source('04_SensisFunctions.R')
source('04_SensisFunctions_Females.R')

sensis_crossed_df<- data.frame()
for (year in as.numeric(colnames(SSE_deathrates_male_df)[-1])){
  sensis_crossed_df <- rbind(sensis_crossed_df,cbind(compute_delta_Ex_crossedsensis(period1 = year-1, period2 = year), year))
}


sensis_crossed_var_df <- (sensis_crossed_df[, c('LE_sensisHump','LE_sensisSenescent','LE_sensisInfant')] - sensis_crossed_df$LE_base)*12
sensis_crossed_var_df$year <- sensis_crossed_df$year
sensis_crossed_var_df$LE_var <- (sensis_crossed_df$LE_after - sensis_crossed_df$LE_base)*12

#crossed sensis plot ----------


sensis_crossed_plot_5mean_m  <- melt(aggregate(sensis_crossed_var_df, list(rep(1:(nrow(sensis_crossed_var_df)%/%5+1),each=5,len=nrow(sensis_crossed_var_df))),mean)[-1], id.vars = 'year' )
sensis_crossed_plot_5mean_m$sex <- 'M'
sensis_crossed_plot_5mean_m$year <- c(seq(1966,2011, 5), 2017)

sensis_crossed_plot_5mean_f  <- melt(aggregate(sensis_crossed_var_df, list(rep(1:(nrow(sensis_crossed_var_df)%/%5+1),each=5,len=nrow(sensis_crossed_var_df))),mean)[-1], id.vars = 'year' )
sensis_crossed_plot_5mean_f$sex <- 'F'
sensis_crossed_plot_5mean_f$year <- seq(1961,2016, 5)

sensis_crossed_plot <- melt(sensis_crossed_var_df, id.vars = 'year' )
sensis_crossed_plot_5mean <- rbind(sensis_crossed_plot_5mean_f, sensis_crossed_plot_5mean_m)

p1 <- ggplot(subset(sensis_crossed_plot,  variable != 'Senescent')) + geom_line(aes(x= year, y= value, group = variable, color= variable)) + #facet_wrap(~ variable) +
  ggtitle('Component Yearly Ex improvement ') + ylab('LE improvement (months)')

p2 <- ggplot(sensis_crossed_plot_5mean) + geom_line(aes(x= year, y= value, group = variable, color= variable, linetype= variable)) +  facet_wrap(~ sex, nrow=2)+
  ggtitle('Component 5Y Ex Improvements')+ ylab('LE improvement (months)') +
  scale_x_continuous(breaks = seq(1961,2016, 5))+
  scale_linetype_manual('LE Improvement',values = c(1,1,1,2), labels = c('Hump', 'Senescent', 'Infant', 'Total'))+
  scale_color_discrete('LE Improvement',  labels = c('Hump', 'Senescent', 'Infant', 'Total'))

png('comp_imp_FR.png', width = 5000, height = 4000, res = 500)
p2
dev.off()

#write.csv(sensis_crossed_df, 'sensis_ex_crossed.csv')
#ggplot(sensis_crossed_df) + geom_line(aes(x= year, y= val)) #+ facet_wrap(~ variable)

#Corrélations ---------------------------------------
png('studyppt.png', height=1500, width=1500, res=200)
plot_selectedyears_allComp(c(2017)) + theme(legend.position = 'none')
dev.off()
plot_qx_crossedsensis(1965,1966,component = 'senescent')

imp_corr <- rbind(melt(cor(dcast(subset(sensis_crossed_plot_5mean, sex=='F'), formula = year ~ variable, value.var = 'value')[,c(2,3,4)])),
                  melt(cor(dcast(subset(sensis_crossed_plot_5mean, sex=='M'), formula = year ~ variable, value.var = 'value')[,c(2,3,4)])))
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


png('corrplot_FR.png', height=1600, width = 3500, res=250)
corr_F
dev.off()


