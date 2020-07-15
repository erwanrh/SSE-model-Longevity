######## SCRIPT MODEL SSE Sensitivities Analysis #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)




####################################################################################-
#-------------------------- SENSITIVITES  ----------------------------
####################################################################################-

splines <- cbind(SSE_males$XX$X1,SSE_males$XX$X2,SSE_males$XX$X3) #Les splines ne sont jamais modifiées
coefs <- matrix(rep(t(SSE_males$coef),times=111),nrow =  111,byrow = T)
eta <- splines*coefs
gamma <- cbind(exp(rowSums(eta[,1:2])), exp(rowSums(eta[,3:27])), exp(rowSums(eta[,28:52]) ))
mHat <- rowSums(gamma)

x11()
plot(mHat,type='l',log='y')
axis(1, at=seq(0,110,5))


# 7.i Modification du coefficient 2  ----------------------------
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

# 7.ii Dataframe des Qx pour chaque valeur de alpha 2  ----------------------------
Qx_Modified_df <- data.frame(matrix(unlist(Qx_Modified_list), ncol=length(Qx_Modified_list), byrow=F))
colnames(Qx_Modified_df) <- names(Qx_Modified_list)
Qx_Modified_df$Age <- 0:110
Qx_Modified_df <- melt(Qx_Modified_df,variable.name = 'Alpha',value.name = 'Qx',id.vars = 'Age',factorsAsStrings = TRUE)
Qx_Modified_df$Alpha <- as.numeric(as.character(Qx_Modified_df$Alpha))
ggplot(Qx_Modified_df) + geom_line( aes(x = Age, y=Qx, group = Alpha, color= Alpha)) +scale_y_continuous(trans='log10') + ggtitle('Courbe de mortalité en fonction de Alpha 2') +
  scale_x_continuous(breaks= seq(0,110,10))


# 7.iii Plot des espérances de vie en fonction du paramètre alpha 2  ----------------------------
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





# 7.iv Sensibilité des coeffs 2 a 2  ----------------------------
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






# 7.v Sensibilite 1 a 1  ----------------------------

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


# 7.vi Dataframe des Qx pour chaque MODIFICATION  ----------------------------
Qx_Modified_df3 <- data.frame(matrix(unlist(Qx_Modified_list3), ncol=length(Qx_Modified_list3), byrow=F))
colnames(Qx_Modified_df3) <- names(Qx_Modified_list3)
Qx_Modified_df3$Age <- 0:110
Qx_Modified_df3 <- melt(Qx_Modified_df3,variable.name = 'ModifiedCoeffs',value.name = 'Qx',id.vars = 'Age',factorsAsStrings = TRUE)

ggplot(Qx_Modified_df3) + geom_line( aes(x = Age, y=Qx, group = ModifiedCoeffs, color= ModifiedCoeffs)) +scale_y_continuous(trans='log10') + ggtitle('Courbe de mortalité en fonction des coefficients modifiés 1 à 1') +
  scale_x_continuous(breaks= seq(0,110,10))+ theme(legend.position="bottom") +  guides(col = guide_legend(nrow = 5))



write.csv2(Study_Coeff_df, 'studycoeff22.csv')
write.csv2(Study_Coeff2_df, 'studycoeff11.csv')

