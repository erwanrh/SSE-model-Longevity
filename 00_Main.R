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
for (pk in pckg$HMDdata){
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
source('02_Plots.R')

####################################################################################-
#---------------------------4. SPLINES  ----------------------------
####################################################################################-

# 4.i Egalite des splines  ----------------------------
SSE_splines_AllGenders <- merge(SSE_splines_female, SSE_splines_male,by = c('age','splinenb','year'), all=T )
SSE_splines_AllGenders$Equal <- SSE_splines_AllGenders$splinevalue.x== SSE_splines_AllGenders$splinevalue.y
if (sum(SSE_splines_AllGenders$Equal) != nrow(SSE_splines_AllGenders)){
  print("Warning : les splines ne sont pas égales pour les Hommes et les Femmes en tout âge et année")
}



# 4.ii Déplacement des splines  ----------------------------

splines_df <- SSE_splines_male[SSE_splines_male$year == 2007,]
splines_df_cast <- dcast(splines_df, age ~ splinenb,value.var = 'splinevalue')
splines_df_cast$age <- NULL

max_index = data.frame()
for ( c in colnames(splines_df_cast)){
  max_index[c,'ArgMax'] <- match(max(splines_df_cast[,c]),splines_df_cast[,c])
  max_index[c,'Max'] <- max(splines_df_cast[,c])
  max_index[c,'SplineNb'] <- as.numeric(c)
}

max_index [max_index$SplineNb %in% 1:2, 'comp']<-as.character(1)
max_index[max_index$SplineNb %in% 3:27, 'comp']<-as.character(2)
max_index[max_index$SplineNb %in% 28:52, 'comp']<-as.character(3)

#Plot de l'âge max et du max par composante et par spline
ggplot(max_index) + geom_point(aes(x=ArgMax, y=Max,color=comp,shape= comp, size=comp)) +
  scale_shape_manual(values=c(16, 16, 16)) + scale_size_manual(values=c(8,5,2)) + 
  scale_color_manual(values = c('#FFD54F','#EA8A75','#7FB3D5'))+
  scale_x_continuous(breaks = seq(0,110,5))


####################################################################################-
#---------------------------5. MUHAT COMPUTATION  ----------------------------
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

#ggsave('plot_allfittedComp.pdf', plot_allFittedcomp, height = 6, width =7 )


#### SENSITIVITY ON COMPONENTS

source('04_SensisFunctions.R')


#Sensitivités -------------------


compute_delta_entropy_sensis(sensis_infant = 0.1, period = 1990)
compute_delta_entropy_sensis(sensis_infant = 0.1, period = 2017)

compute_delta_Ex_sensis(sensis_infant =  0.1, sensis_hump = 0.1, sensis_senescent = 0.1)


plot_lx_sensis(sensis_infant =  100, sensis_hump = 0.1, sensis_senescent = 0)
plot_qx_sensis(sensis_infant =  100, period = 1965)

sensis_ex_df <- data.frame()
sensis_entropy_df <- data.frame()

for(delta in c(0.1, 0.2, 0.5, 0.7, 1, 5,10, 100)){
  for (year in 1960:2017){
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'Infant'] <- compute_delta_Ex_sensis(sensis_infant = delta, period = year)
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'Hump'] <- compute_delta_Ex_sensis(sensis_hump = delta, period = year)
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'Senescent'] <- compute_delta_Ex_sensis(sensis_senescent = delta, period = year)
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'variation'] <- delta
    sensis_ex_df[as.character(paste(year,delta, sep = '_')), 'year'] <- year
  }
  
}

write.csv2(sensis_ex_df, 'sensis_life_expectancy2.csv')
write.csv2(sensis_entropy_df, 'sensis_entropy.csv')
#ggsave('plot_allfittedComp_Sensis.pdf', plot_allFittedcomp, height = 6, width =7 )



# Crossed sensitivities ---------------------------------------------------

sensis_crossed_df<- data.frame()

for (year in as.numeric(colnames(SSE_deathrates_male_df)[-1])){
  sensis_crossed_df <- rbind(sensis_crossed_df,compute_delta_Ex_crossedsensis(period1 = year-1, period2 = year))
}
sensis_crossed_df <- sensis_crossed_df*100
sensis_crossed_df$year <- as.numeric(colnames(SSE_deathrates_male_df)[-1])
colnames(sensis_crossed_df) <- c('Hump', 'Senescent', 'Infant','LE_before','LE_after', 'year')
colnames(sensis_crossed_df) <- c('val', 'year')

write.csv2(sensis_crossed_df, 'sensis_ex_crossed.csv')

sensis_crossed_plot <- melt(sensis_crossed_df, id.vars = 'year' )

ggplot(subset(sensis_crossed_plot,  variable != 'Senescent')) + geom_line(aes(x= year, y= value, group = variable, color= variable)) #+ facet_wrap(~ variable)
ggplot(sensis_crossed_df) + geom_line(aes(x= year, y= val)) #+ facet_wrap(~ variable)

mean(sensis_crossed_df$val)

plot_selectedyears_allComp(c(1965,1966))
