#install.packages('HMDHFDplus')
rm(list=ls())

library(HMDHFDplus)
library(Matrix)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
      
library("zoo")         

setwd('/Users/erwanrahis/Desktop//SSE_tests')


#MORTSMOOTH 
pckg<-read.csv('pckSSE.csv')
for (pk in pckg$HMDdata){
  source(paste0('MortalitySmooth-master/MortalitySmooth-master/R/',pk,'.R'))  
}

#MORTHUMP
source('SSE_fit.R')




## FIT DU MODELE SUR PLUSIEURS ANNEES 
country_code <- 'FR'
year_min <- 1960
year_max <- 2005

SSE_coeffcients_males <-list()
SSE_coeffcients_females <- list()
SSE_deathrates_male <- list()
SSE_deathrates_female <- list()
SSE_splines_male <- NULL
SSE_splines_female <-NULL



#FIT du modèle sur chaque année entre 1900 et 2017
for (year in seq(year_min,year_max,5)){
  Data_males <- HMD2MH(country=country_code,year=year, sex='males',path='Data_France',xtra=TRUE)
  Data_females <- HMD2MH(country=country_code,year=year, sex='females',path='Data_France',xtra=TRUE)
  # Expo
  Data_males$n[ Data_males$n==0] <-0.01
  Data_females$n[ Data_females$n==0] <-0.01
  #Deaths
  Data_females$d <- as.integer(Data_females$d)
  Data_males$d <- as.integer(Data_males$d)
  #Rates
  Data_females$m[ Data_females$d ==0] <-  0
  Data_males$m[ Data_males$d ==0] <-  0
  
  
  if (!any(is.na(Data_males)) & !any(is.na(Data_females)) ){
    #Fit
    SSE_males<- morthump(data=Data_males, model='sse')
    SSE_females<- morthump(data=Data_females, model='sse')
    
    
    #Récupération des coefficients alpha
    SSE_coeffcients_males[as.character(year)] <-as.data.frame(as.matrix(coef(SSE_males)))
    SSE_coeffcients_females[as.character(year)] <- as.data.frame(as.matrix(coef(SSE_females)))
    
    #Récupération des taux fittés
    SSE_deathrates_male[as.character(year)] <- as.data.frame(SSE_males$mhat$mhat3[,1]+SSE_males$mhat$mhat1[,1]+SSE_males$mhat$mhat2[,1])
    SSE_deathrates_female[as.character(year)] <- as.data.frame(SSE_females$mhat$mhat3[,1]+SSE_females$mhat$mhat1[,1]+SSE_females$mhat$mhat2[,1])
    
    #Récupération des composantes
    temp_SSE_splines_male <- melt(cbind(SSE_males$XX$X1,SSE_males$XX$X2,SSE_males$XX$X3), value.name = 'splinevalue',varnames = c('age','splinenb'))
    temp_SSE_splines_male$year <- year
    
    temp_SSE_splines_female <- melt(cbind(SSE_females$XX$X1,SSE_females$XX$X2,SSE_females$XX$X3), value.name = 'splinevalue',varnames = c('age','splinenb'))
    temp_SSE_splines_female$year <- year
    
    SSE_splines_male<- rbind(SSE_splines_male,temp_SSE_splines_male)
    SSE_splines_female<- rbind(SSE_splines_female,temp_SSE_splines_female)
    }
  
}



#Plot des courbes de mortalité
SSE_deathrates_male_df <-  t(data.frame(matrix(unlist(SSE_deathrates_male), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_deathrates_female_df <-  t(data.frame(matrix(unlist(SSE_deathrates_female), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))

DeathRates_plot_males <- melt(SSE_deathrates_male_df,varnames = c('Age','Year'),value.name = 'coeff',)
DeathRates_plot_females <- melt(SSE_deathrates_female_df,varnames = c('Age','Year'),value.name = 'coeff',)
DeathRates_plot_males$gender <- 'M'
DeathRates_plot_females$gender <- 'F'
DeathRates_plot <- rbind(DeathRates_plot_females,DeathRates_plot_males)
DeathRates_plot$Age <- as.numeric(DeathRates_plot$Age)

DeathRatesplot <- ggplot(DeathRates_plot)+geom_line(aes(x=Age, y=coeff, group = Year,color=Year))+ ggtitle(paste0('SSE Mortality Rates - ',country_code)) +facet_wrap(gender~.,scales = 'free')+
  scale_y_continuous(trans='log2') + scale_x_continuous(breaks = seq(0,110,5), labels= seq(0,110,5) )
ggsave(paste0('drplot_',country_code,'.pdf'),DeathRatesplot,width = 20, height = 12)


#Export des taux de mortalité en CSV
write.csv(SSE_deathrates_female_df, paste0('qx_F_', country_code,'.csv'))
write.csv(SSE_deathrates_male_df, paste0('qx_M_', country_code,'.csv'))



#Mise en forme des données des coefficients pour plot
SSE_coeffcients_males_df <- t(data.frame(matrix(unlist(SSE_coeffcients_males), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_coeffcients_females_df <- t(data.frame(matrix(unlist(SSE_coeffcients_females), nrow=length(SSE_coeffcients_females), byrow=T),row.names = names(SSE_coeffcients_females)))

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
ggsave(paste0('coefplots_',country_code,'.pdf'),coefplots,width = 20, height = 12)







#GRAPH OF MUHAT OUTPUTS
png(paste0('TestSSE_outputs_',country_code,'.png'), units="in", width=10, height=9, res=300 )

plot(SSE_males$mhat$mhat3[,1]+SSE_males$mhat$mhat1[,1]+SSE_males$mhat$mhat2[,1],col='red',log='y',type='l',ylim = c(0.00001,1),ylab = 'mhat',xlab='age')
lines(SSE_males$mhat$mhat1[,1],type='l',col='blue',lty=2)
lines(SSE_males$mhat$mhat2[,1],col='skyblue',lty=2)
lines(SSE_males$mhat$mhat3[,1],col='darkblue',lty=2)
points(SSE_males$data$m,col='pink')
legend(0,1,legend=c('Composante 1','Composante 2','Composante 3','1 + 2 + 3','Data'),col = c('blue','skyblue','darkblue','red','pink'),lty=c(2,2,2,1,1))
title(paste0('Components of SSE model - ',as.character(year),' - ', country_code))

dev.off()




#Test d'égalité des splines homme et femme et par années
SSE_splines_AllGenders <- merge(SSE_splines_female, SSE_splines_male,by = c('age','splinenb','year'), all=T )
SSE_splines_AllGenders$Equal <- SSE_splines_AllGenders$splinevalue.x== SSE_splines_AllGenders$splinevalue.y
if (sum(SSE_splines_AllGenders$Equal) != nrow(SSE_splines_AllGenders)){
  print("Warning : les splines ne sont pas égales pour les Hommes et les Femmes en tout âge et année")
}

#Graph des élements des B-splines
splines_plot <- SSE_splines_male
splines_plot [splines_plot$splinenb %in% 1:2, 'comp']<-'c1'
splines_plot[splines_plot$splinenb %in% 3:27, 'comp']<-'c2'
splines_plot[splines_plot$splinenb %in% 28:52, 'comp']<-'c3'


splines_plot_2 <- dcast(splines_plot,age+splinenb~year,value.var = 'splinevalue') 

splines_ggplot <- ggplot(splines_plot,aes(x=age, y=splinevalue, color=year, group=year))+ geom_line() + facet_wrap(.~splinenb+ comp) + ggtitle(paste0('Splines - ', country_code))
ggsave(paste0('spline_comp_',country_code,'.pdf'),plot = splines_ggplot,width = 15, height = 15)
   

# Graph de tous les coefficients pour comparaison avec les splines
coefs_plots_52gridM <- ggplot(Coeff_plot_males,aes(x=Year, y=coeff, color=year, group=year))+ geom_line() + facet_wrap(.~Coeff.num) + ggtitle (paste0('Coefficients for splines Males - ', country_code))
ggsave(paste0('coeffs_52grid_M',country_code,'.pdf'),plot = coefs_plots_52gridM,width = 15, height = 15) 

coefs_plots_52gridF <- ggplot(Coeff_plot_females,aes(x=Year, y=coeff, color=year, group=year))+ geom_line() + facet_wrap(.~Coeff.num) + ggtitle (paste0('Coefficients for splines Females - ', country_code))
ggsave(paste0('coeffs_52grid_F',country_code,'.pdf'),plot = coefs_plots_52gridF,width = 15, height = 15) 









######  LIFE EXPECTANCIES
source('LifeExpectancies.R')
QX_Ex_Male <- dcast(DeathRates_plot_males,formula = Age ~ Year,value.var = 'coeff')
QX_Ex_Male$Age <- NULL

LE_FR_df <- data.frame()

for (i in c(1960,1965,1970,1975,1980,1985,1990,1995,2000,2005)){
  LE_FR_df['LE_birth_base',as.character(i)] <- LE_period_Model(0,i,QX_Ex_Male)
}

write.csv2(LE_FR_df, 'lefr.csv')


#Sensibilité des coefficients
splines <- cbind(SSE_males$XX$X1,SSE_males$XX$X2,SSE_males$XX$X3) #Les splines ne sont jamais modifiées
coefs <- matrix(rep(t(SSE_males$coef),times=111),nrow =  111,byrow = T)
eta <- splines*coefs
gamma <- cbind(exp(rowSums(eta[,1:2])), exp(rowSums(eta[,3:27])), exp(rowSums(eta[,28:52]) ))
mHat <- rowSums(gamma)

x11()
plot(mHat,type='l',log='y')
title(paste0('Courbe de mortalité - ', year))
axis(1, at=seq(0,110,5))

#########Modification du coefficient 2
coeff2_lifeExp <- data.frame()
Qx_Modified_list <- {}

for (x in seq(0,6,1)){
  coef1 <- SSE_males$coef
  coef1[2] <-  x
  coeffs_modification1 <- matrix(rep(t(coef1),times=111),nrow =  111,byrow = T)
  
  eta_modification1 <- splines*coeffs_modification1
  gamma_modification1 <- cbind(exp(rowSums(eta_modification1[,1:2])), exp(rowSums(eta_modification1[,3:27])), exp(rowSums(eta_modification1[,28:52]) ))
  mHat_modification1 <- rowSums(gamma_modification1)
  Qx_Modified<- data.frame('2007'=mHat_modification1,check.names = F)
  Qx_Modified_list[as.character(x)] <-  Qx_Modified
  
  coeff2_lifeExp[as.character(x),'LifeEx']<-LE_period_Model(0,2007,Qx_Modified)
  coeff2_lifeExp[as.character(x),'Alpha2']<- x
}

## Dataframe des Qx pour chaque valeur de alpha 2
Qx_Modified_df <- data.frame(matrix(unlist(Qx_Modified_list), ncol=length(Qx_Modified_list), byrow=F))
colnames(Qx_Modified_df) <- names(Qx_Modified_list)
Qx_Modified_df$Age <- 0:110
Qx_Modified_df <- melt(Qx_Modified_df,variable.name = 'Alpha',value.name = 'Qx',id.vars = 'Age',factorsAsStrings = TRUE)
Qx_Modified_df$Alpha <- as.numeric(as.character(Qx_Modified_df$Alpha))
ggplot(Qx_Modified_df) + geom_line( aes(x = Age, y=Qx, group = Alpha, color= Alpha)) +scale_y_continuous(trans='log10') + ggtitle('Courbe de mortalité en fonction de Alpha 2') +
  scale_x_continuous(breaks= seq(0,110,10))





#Plot des espérances de vie en fonction du paramètre alpha 2
plot(coeff2_lifeExp$Alpha2, coeff2_lifeExp$LifeEx, xlab = 'Alpha2', ylab='Life Expectancy')

#Espérance de vie avant modificiations
LE_period_Model(0,2007,QX_Ex_Male)

#x11()  
plot(mHat_modification1,type='l',log='y', ylab='Log-Mortality',xlab='Age')



#### MODIFICATION DES PARAMETRES 29/30

coef1 <- SSE_males$coef
coef1[44] <-  0
coeffs_modification1 <- matrix(rep(t(coef1),times=111),nrow =  111,byrow = T)

eta_modification1 <- splines*coeffs_modification1
gamma_modification1 <- cbind(exp(rowSums(eta_modification1[,1:2])), exp(rowSums(eta_modification1[,3:27])), exp(rowSums(eta_modification1[,28:52]) ))
mHat_modification1 <- rowSums(gamma_modification1)
Qx_Modified<- data.frame('1990'=mHat_modification1,check.names = F)

LE_period_Model(0,1990,QX_Ex_Male)
LE_period_Model(0,1990,Qx_Modified)

plot(mHat_modification1,type='l',log='y', ylab='Log-Mortality',xlab='Age', col='green')
lines(mHat)
axis(1, at = seq(0,110,10))





