######## SCRIPT MODELE SSE ####### #
# Camarda, 2016, Sum Of Smooth Exponentials

####1. INITIALISATION####
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



####2. MODEL FITTING ####
country_code <- 'CHL'
year_min <- 2005
year_max <- 2008
years <- c(1992,1994,1995,1996,1998,1999,2000,2001,2002,2003,2005,2006,2007)
   
SSE_coeffcients_males <-list()
SSE_coeffcients_females <- list()
SSE_deathrates_male <- list()
SSE_deathrates_female <- list()
SSE_splines_male <- NULL
SSE_splines_female <-NULL


for (year in years){
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

# DATA FRAMES Coefficients et Death Rates
SSE_coeffcients_males_df <- t(data.frame(matrix(unlist(SSE_coeffcients_males), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_coeffcients_females_df <- t(data.frame(matrix(unlist(SSE_coeffcients_females), nrow=length(SSE_coeffcients_females), byrow=T),row.names = names(SSE_coeffcients_females)))
SSE_deathrates_male_df <-  t(data.frame(matrix(unlist(SSE_deathrates_male), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))
SSE_deathrates_female_df <-  t(data.frame(matrix(unlist(SSE_deathrates_female), nrow=length(SSE_coeffcients_males), byrow=T),row.names = names(SSE_coeffcients_males)))



# 3. PLOTS ####

# 3.i Plot des courbes de mortalité ####
DeathRates_plot_males <- melt(SSE_deathrates_male_df,varnames = c('Age','Year'),value.name = 'coeff',)
DeathRates_plot_females <- melt(SSE_deathrates_female_df,varnames = c('Age','Year'),value.name = 'coeff',)
DeathRates_plot_males$gender <- 'M'
DeathRates_plot_females$gender <- 'F'
DeathRates_plot <- rbind(DeathRates_plot_females,DeathRates_plot_males)
DeathRates_plot$Age <- as.numeric(DeathRates_plot$Age)

DeathRatesplot <- ggplot(DeathRates_plot)+geom_line(aes(x=Age, y=coeff, group = Year,color=Year))+ ggtitle(paste0('SSE Mortality Rates - ',country_code)) +facet_wrap(gender~.,scales = 'free')+
  scale_y_continuous(trans='log2') + scale_x_continuous(breaks = seq(0,110,5), labels= seq(0,110,5) )

ggsave(paste0('drplot_',country_code,'.pdf'),DeathRatesplot,width = 20, height = 12)


# 3.ii Plot des coefficients ####
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




#3.iii Plot des muhat ####
png(paste0('TestSSE_outputs_',country_code,'.png'), units="in", width=10, height=9, res=300 )

plot(SSE_males$mhat$mhat3[,1]+SSE_males$mhat$mhat1[,1]+SSE_males$mhat$mhat2[,1],col='red',log='y',type='l',ylim = c(0.00001,1),ylab = 'mhat',xlab='age')
lines(SSE_males$mhat$mhat1[,1],type='l',col='blue',lty=2)
lines(SSE_males$mhat$mhat2[,1],col='skyblue',lty=2)
lines(SSE_males$mhat$mhat3[,1],col='darkblue',lty=2)
points(SSE_males$data$m,col='pink')
legend(0,1,legend=c('Composante 1','Composante 2','Composante 3','1 + 2 + 3','Data'),col = c('blue','skyblue','darkblue','red','pink'),lty=c(2,2,2,1,1))
title(paste0('Components of SSE model - ',as.character(year),' - ', country_code))

dev.off()


# 3.iv Graph des élements des B-splines ####
splines_plot <- SSE_splines_male
splines_plot [splines_plot$splinenb %in% 1:2, 'comp']<-'c1'
splines_plot[splines_plot$splinenb %in% 3:27, 'comp']<-'c2'
splines_plot[splines_plot$splinenb %in% 28:52, 'comp']<-'c3'

splines_ggplot <- ggplot(splines_plot,aes(x=age, y=splinevalue, color=year, group=year))+ geom_line() + facet_wrap(.~splinenb+ comp) + ggtitle(paste0('Splines - ', country_code))
ggsave(paste0('spline_comp_',country_code,'.pdf'),plot = splines_ggplot,width = 15, height = 15)


# 3.v Graph de tous les coefficients pour comparaison avec les splines ####
coefs_plots_52gridM <- ggplot(Coeff_plot_males,aes(x=Year, y=coeff, color=year, group=year))+ geom_line() + facet_wrap(.~Coeff.num) + ggtitle (paste0('Coefficients for splines Males - ', country_code))
ggsave(paste0('coeffs_52grid_M',country_code,'.pdf'),plot = coefs_plots_52gridM,width = 15, height = 15) 

coefs_plots_52gridF <- ggplot(Coeff_plot_females,aes(x=Year, y=coeff, color=year, group=year))+ geom_line() + facet_wrap(.~Coeff.num) + ggtitle (paste0('Coefficients for splines Females - ', country_code))
ggsave(paste0('coeffs_52grid_F',country_code,'.pdf'),plot = coefs_plots_52gridF,width = 15, height = 15) 



# 4.  SPLINES ####

# 4.i Egalite des splines ####
SSE_splines_AllGenders <- merge(SSE_splines_female, SSE_splines_male,by = c('age','splinenb','year'), all=T )
SSE_splines_AllGenders$Equal <- SSE_splines_AllGenders$splinevalue.x== SSE_splines_AllGenders$splinevalue.y
if (sum(SSE_splines_AllGenders$Equal) != nrow(SSE_splines_AllGenders)){
  print("Warning : les splines ne sont pas égales pour les Hommes et les Femmes en tout âge et année")
}



# 4.ii Déplacement des splines ####

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



# 5. CALCUL DES MUHAT #####
plot(SSE_males$data$m,col='pink',type='o',log='y',ylim = c(0.00001,1))

lines(exp(SSE_males$XX$X1%*%(coef(SSE_males)[1:2])),col='red')
lines(exp(SSE_males$XX$X2%*%(coef(SSE_males)[3:27])),col='green')
lines(exp(SSE_males$XX$X3%*%(coef(SSE_males)[28:52])),col='blue')

hump <- melt(exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
hump$Comp <- '3'
ggplot(hump) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

infant <- melt(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2,])),varnames = c('Age','Year'),value.name = 'm')
infant$Comp <- '1'
ggplot(infant) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

senescent <- melt(exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
senescent$Comp <- '2'
ggplot(senescent) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

MU <- hump
MU$m <- hump$m+infant$m+senescent$m
MU$Comp <-'All'
ggplot(MU) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

ALL_plot_df <- rbind(hump, infant, senescent, MU)
ALL_plot_df$Year <-as.character(ALL_plot_df$Year)

ggplot(subset(ALL_plot_df, Year %in% c(2000,2004,1996))) + 
  geom_line( aes(x=Age, y=m, group= interaction(Year,Comp), color=Comp,linetype=Year))+
  scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ')  +ggtitle('SSE Composantes')




# 6.  LIFE EXPECTANCIES ####
source('LifeExpectancies.R')
QX_Ex_Male <- dcast(DeathRates_plot_males,formula = Age ~ Year,value.var = 'coeff')
QX_Ex_Male$Age <- NULL

LE_CHL_df <- data.frame()

for (i in c(1992,1996,2000,2005,2007)){
  LE_CHL_df['LE_birth_base',as.character(i)] <- LE_period_Model(0,i,QX_Ex_Male)
}

write.csv(LE_CHL_df,'lechl.csv')



# 7. SENSIBILITES ####
splines <- cbind(SSE_males$XX$X1,SSE_males$XX$X2,SSE_males$XX$X3) #Les splines ne sont jamais modifiées
coefs <- matrix(rep(t(SSE_males$coef),times=111),nrow =  111,byrow = T)
eta <- splines*coefs
gamma <- cbind(exp(rowSums(eta[,1:2])), exp(rowSums(eta[,3:27])), exp(rowSums(eta[,28:52]) ))
mHat <- rowSums(gamma)

x11()
plot(mHat,type='l',log='y')
axis(1, at=seq(0,110,5))


# 7.i Modification du coefficient 2 ####
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

# 7.ii Dataframe des Qx pour chaque valeur de alpha 2 ####
Qx_Modified_df <- data.frame(matrix(unlist(Qx_Modified_list), ncol=length(Qx_Modified_list), byrow=F))
colnames(Qx_Modified_df) <- names(Qx_Modified_list)
Qx_Modified_df$Age <- 0:110
Qx_Modified_df <- melt(Qx_Modified_df,variable.name = 'Alpha',value.name = 'Qx',id.vars = 'Age',factorsAsStrings = TRUE)
Qx_Modified_df$Alpha <- as.numeric(as.character(Qx_Modified_df$Alpha))
ggplot(Qx_Modified_df) + geom_line( aes(x = Age, y=Qx, group = Alpha, color= Alpha)) +scale_y_continuous(trans='log10') + ggtitle('Courbe de mortalité en fonction de Alpha 2') +
  scale_x_continuous(breaks= seq(0,110,10))


# 7.iii Plot des espérances de vie en fonction du paramètre alpha 2 ####
plot(coeff2_lifeExp$Alpha2, coeff2_lifeExp$LifeEx, xlab = 'Alpha2', ylab='Life Expectancy')

#Espérance de vie avant modificiations
LE_period_Model(0,2007,QX_Ex_Male)

#x11()
plot(mHat_modification1,type='l',log='y', ylab='Log-Mortality',xlab='Age')



#### MODIFICATION DES PARAMETRES 29/30

coef1 <- SSE_males$coef
coef1[28:52] <-  0
coeffs_modification1 <- matrix(rep(t(coef1),times=111),nrow =  111,byrow = T)

eta_modification1 <- splines*coeffs_modification1
gamma_modification1 <- cbind(exp(rowSums(eta_modification1[,1:2])), exp(rowSums(eta_modification1[,3:27])), exp(rowSums(eta_modification1[,28:52]) ))
mHat_modification1 <- rowSums(gamma_modification1)
Qx_Modified<- data.frame('2007'=mHat_modification1,check.names = F)

LE_period_Model(0,2007,QX_Ex_Male)
LE_period_Model(0,2007,Qx_Modified)

plot(mHat_modification1,type='l',log='y', ylab='Log-Mortality',xlab='Age', col='green')
lines(mHat)
axis(1, at = seq(0,110,10))



#Plot de l'exposition
plot(SSE_males$data$n,ylab = 'Exposition', xlab='Age', type='l')
title('Exposition en fonction de l‘âge')





# 7.iv Sensibilité des coeffs 2 a 2 ####
Coefficient_variation <- data.frame()
for (coeff in row.names(SSE_coeffcients_males_df)){
  
  Coefficient_variation[coeff,'CoeffVar'] <- abs(sd(SSE_coeffcients_males_df[coeff,])/ mean(SSE_coeffcients_males_df[coeff,]))
  Coefficient_variation[coeff,'Mean'] <- mean(SSE_coeffcients_males_df[coeff,])
  Coefficient_variation[coeff,'Std'] <- sd(SSE_coeffcients_males_df[coeff,])
  
}


Study_Coeff_df <- data.frame()
Qx_Modified_list2 <- {}

for (i2 in 3:27){
  i3 <- i2+25
  
  coef1 <- SSE_males$coef
  coef1[i2] <-  coef1[i2]*(1-Coefficient_variation[i2,1])
  coef1[i3]<-  coef1[i2]*(1-Coefficient_variation[i3,1])
  
  coeffs_modification1 <- matrix(rep(t(coef1),times=111),nrow =  111,byrow = T)
  
  eta_modification1 <- splines*coeffs_modification1
  gamma_modification1 <- cbind(exp(rowSums(eta_modification1[,1:2])), exp(rowSums(eta_modification1[,3:27])), exp(rowSums(eta_modification1[,28:52]) ))
  mHat_modification1 <- rowSums(gamma_modification1)
  Qx_Modified<- data.frame('2007'=mHat_modification1,check.names = F)
  Qx_Modified_list2[paste0(as.character(i2),('-'),as.character(i3))] <-  Qx_Modified
  
  Study_Coeff_df[paste0(as.character(i2),('-'),as.character(i3)),'IndiceComposante2']<- i2
  Study_Coeff_df[paste0(as.character(i2),('-'),as.character(i3)),'IndiceComposante3']<- i3
  Study_Coeff_df[paste0(as.character(i2),('-'),as.character(i3)),'LifeEx']<-LE_period_Model(0,2007,Qx_Modified)
  Study_Coeff_df[paste0(as.character(i2),('-'),as.character(i3)),'VariationEx']<- LE_period_Model(0,2007,QX_Ex_Male)-LE_period_Model(0,2007,Qx_Modified)
  
}

## Dataframe des Qx pour chaque modification
Qx_Modified_df2 <- data.frame(matrix(unlist(Qx_Modified_list2), ncol=length(Qx_Modified_list2), byrow=F))
colnames(Qx_Modified_df2) <- names(Qx_Modified_list2)
Qx_Modified_df2$Age <- 0:110
Qx_Modified_df2 <- melt(Qx_Modified_df2,variable.name = 'ModifiedCoeffs',value.name = 'Qx',id.vars = 'Age',factorsAsStrings = TRUE)

ggplot(Qx_Modified_df2) + geom_line( aes(x = Age, y=Qx, group = ModifiedCoeffs, color= ModifiedCoeffs)) +scale_y_continuous(trans='log10') + ggtitle('Courbe de mortalité en fonction des coefficients modifiés 2 à 2') +
  scale_x_continuous(breaks= seq(0,110,10))+ theme(legend.position="bottom")+  guides(col = guide_legend(nrow = 5))






# 7.v Sensibilite 1 a 1 ####

Study_Coeff2_df <- data.frame()
Qx_Modified_list3 <- {}

for (i2 in 1:52){

  coef1 <- SSE_males$coef
  coef1[i2] <-  coef1[i2]*(1-Coefficient_variation[i2,1])

  
  coeffs_modification1 <- matrix(rep(t(coef1),times=111),nrow =  111,byrow = T)
  
  eta_modification1 <- splines*coeffs_modification1
  gamma_modification1 <- cbind(exp(rowSums(eta_modification1[,1:2])), exp(rowSums(eta_modification1[,3:27])), exp(rowSums(eta_modification1[,28:52]) ))
  mHat_modification1 <- rowSums(gamma_modification1)
  Qx_Modified<- data.frame('2007'=mHat_modification1,check.names = F)
  Qx_Modified_list3[as.character(i2)] <-  Qx_Modified
  
  Study_Coeff2_df[as.character(i2),'IndiceComposante']<- i2
  Study_Coeff2_df[as.character(i2),'LifeEx']<-LE_period_Model(0,2007,Qx_Modified)
  Study_Coeff2_df[as.character(i2),'VariationEx']<- LE_period_Model(0,2007,QX_Ex_Male)-LE_period_Model(0,2007,Qx_Modified)
  
}


# 7.vi Dataframe des Qx pour chaque MODIFICATION ####
Qx_Modified_df3 <- data.frame(matrix(unlist(Qx_Modified_list3), ncol=length(Qx_Modified_list3), byrow=F))
colnames(Qx_Modified_df3) <- names(Qx_Modified_list3)
Qx_Modified_df3$Age <- 0:110
Qx_Modified_df3 <- melt(Qx_Modified_df3,variable.name = 'ModifiedCoeffs',value.name = 'Qx',id.vars = 'Age',factorsAsStrings = TRUE)

ggplot(Qx_Modified_df3) + geom_line( aes(x = Age, y=Qx, group = ModifiedCoeffs, color= ModifiedCoeffs)) +scale_y_continuous(trans='log10') + ggtitle('Courbe de mortalité en fonction des coefficients modifiés 1 à 1') +
  scale_x_continuous(breaks= seq(0,110,10))+ theme(legend.position="bottom") +  guides(col = guide_legend(nrow = 5))



write.csv2(Study_Coeff_df, 'studycoeff22.csv')
write.csv2(Study_Coeff2_df, 'studycoeff11.csv')




#8. SENSIBILITE DES COMPOSANTES ####
#Admettons que nous gardons la bosse des accidents constantes car nous n'avons aucune idée sur son évolution. 
#Nous faisons évoluer la mortalité infantile.

plot(SSE_males$data$m,col='pink',type='o',log='y',ylim = c(0.00001,1))
lines(exp(SSE_males$XX$X1%*%(coef(SSE_males)[1:2])),col='red')
lines(exp(SSE_males$XX$X2%*%(coef(SSE_males)[3:27])),col='green')
lines(exp(SSE_males$XX$X3%*%(coef(SSE_males)[28:52])),col='blue')

hump2 <- melt(exp(SSE_males$XX$X3%*%(SSE_coeffcients_males_df[28:52,])),varnames = c('Age','Year'),value.name = 'm')
hump2$Comp <- '3'
ggplot(hump) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

infant <- melt(exp(SSE_males$XX$X1%*%(SSE_coeffcients_males_df[1:2,])),varnames = c('Age','Year'),value.name = 'm')
infant$Comp <- '1'
ggplot(infant) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

senescent <- melt(exp(SSE_males$XX$X2%*%(SSE_coeffcients_males_df[3:27,])),varnames = c('Age','Year'),value.name = 'm')
senescent$Comp <- '2'
ggplot(senescent) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

MU <- hump
MU$m <- hump$m+infant$m+senescent$m
MU$Comp <-'All'
ggplot(MU) + geom_line( aes(x=Age, y=m, group=Year, color=Year))+scale_y_continuous(trans='log10') 

ALL_plot_df <- rbind(hump, infant, senescent, MU)
ALL_plot_df$Year <-as.character(ALL_plot_df$Year)

ggplot(subset(ALL_plot_df, Year %in% c(1996,2000))) + 
  geom_line( aes(x=Age, y=m, group= interaction(Year,Comp), color=Comp,linetype=Year))+
  scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ')  +ggtitle('SSE Composantes')



