######## SCRIPT MODEL SSE data for France #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)


####################################################################################-
#--------------------------- ROBUSTESS 1 ----------------------------
####################################################################################-

#Fit des données de 2004 avec évolution de la composante Hump en 2005
robustness_20042005_list <- compute_fitted_crossedsensis_Hump(2004, 2005)

robustness_2004_df <- robustness_20042005_list[['df_base']]
robustness_20042005_df <- robustness_20042005_list[['df_sensis']]

#Espérance de vie en 2004 et en 2004 et avec la composante accident de 2005
LE_period_Model(Age = 0, Period = 'FittedCurve', QxModel = robustness_20042005_df )
LE_period_Model(Age = 0, Period = 'FittedCurve', QxModel = robustness_2004_df )

#Graphique du fit 2004 avec la composante accident 2005
png('crossedhumpMAR2004.png', width= 3500, height=2000, res=400)
plot_qx_crossedsensis(component = 'hump', period1 = 2004, period2 =  2005)
dev.off()

#Reconstruction des données brutes 
#Données brutes de 2004
Robustness_data <- subset(All_interpolated_data, sex== 'M' & year == 2004, select = c('Age', 'qx_male', 'qx_malenx', 'qx_maleqx') )
colnames(Robustness_data) <- c('x','d','n','m')
Robustness_data <- Robustness_data[order(Robustness_data$x),]
#Remplacement des mx par les mx 2004 avec comp accident 2005
Robustness_data$m <- robustness_20042005_df$FittedCurve
#Calcul du nombre de décès avec d = m * expo -> Exposition constante
Robustness_data$d <- Robustness_data$m * Robustness_data$n

#Fit du modèle sur les nouvelles données 
SSE_robustness<- morthump(data=Robustness_data, model='sse')
#Récupération des coefficients alpha
SSE_coeffcients_robustness <-as.data.frame(as.matrix(coef(SSE_robustness)))
#Récupération des taux fittés
SSE_deathrates_robustness <- as.data.frame(SSE_robustness$mhat$mhat3[,1]+SSE_robustness$mhat$mhat1[,1]+SSE_robustness$mhat$mhat2[,1])
row.names(SSE_deathrates_robustness) <- 0:85

#Plot des composantes et de la courbé fittée
hump_r <- melt(exp(SSE_robustness$XX$X3%*%(SSE_coeffcients_robustness[23:42,])),varnames = c('Age','Year'),value.name = 'm')
hump_r$Comp <- 'Hump'
senescent_r <- melt(exp(SSE_robustness$XX$X2%*%(SSE_coeffcients_robustness[3:22,])),varnames = c('Age','Year'),value.name = 'm')
senescent_r$Comp <- 'Senescent'
infant_r <-  melt(exp(SSE_robustness$XX$X1%*%(SSE_coeffcients_robustness[1:2,])), varnames = c('Age','Year'),value.name = 'm')
infant_r$Comp <- 'Infant'
#PLot du comportement des Qx Fittés
MU_r2 <- hump_r
MU_r2$m <- hump_r$m+infant_r$m+senescent_r$m
MU_r2$Comp <-'FittedCurve'

#Plot de toutes les composantes + des qx fittés
ALL_plot_robust <- rbind(infant_r, hump_r, senescent_r, MU_r2)
ALL_plot_robust$Year <-'2004_newfit'
ALL_plot_robust$Age <- 0:85

ALL_plot_robust2 <- melt(robustness_20042005_df[,-6], id.vars = c('Age'), value.name = 'm', varnames = c('Infant','Senescent', 'Hump', 'FittedCurve'), variable.name = 'Comp')
ALL_plot_robust2$Year <- '2004_before'

ALL_plot_robust12 <- rbind(ALL_plot_robust, ALL_plot_robust2)

plot_robust <- ggplot(subset(ALL_plot_robust12)) + 
  geom_line( aes(x=Age, y=m, group= interaction(Year,Comp), color=Comp,linetype=Year))+
  scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ')  +ggtitle('SSE Composantes')

#EV du nouveau fit 
LE_period_Model(Age = 0, Period = 'm', QxModel = MU_r2 )


#Graphique du fit 2004 avec la composante accident 2005
png('robustMAR.png', width= 3500, height=2000, res=400)
plot_robust
dev.off()




####################################################################################-
#--------------------------- ROBUSTESS 2 ----------------------------
####################################################################################-

#Fit des données de 2016
robustness_2016_list <- compute_fitted_sensis(period = 2016, sensis_hump = 0.05 )

robustness_2016_df <- robustness_2016_list[['df_base']]
robustness_2016sensis_df <- robustness_2016_list[['df_sensis']]

#Espérance de vie en 2004 et en 2004 et avec la composante accident de 2005
LE_period_Model(Age = 0, Period = 'FittedCurve', QxModel = robustness_2016_df  )
LE_period_Model(Age = 0, Period = 'FittedCurve', QxModel = robustness_2016sensis_df )

#Graphique du fit 2004 avec la composante accident 2005
png('sensisseneMAR2016.png', width= 3500, height=2000, res=400)
plot_qx_sensis(sensis_hump = 0.05, period = 2016)
dev.off()

#Reconstruction des données brutes 
#Données brutes de 2004
Robustness2_data <- subset(All_interpolated_data, sex== 'M' & year == 2016, select = c('Age', 'qx_male', 'qx_malenx', 'qx_maleqx') )
colnames(Robustness2_data) <- c('x','d','n','m')
Robustness2_data <- Robustness2_data[order(Robustness2_data$x),]
#Remplacement des mx par les mx 2004 avec comp accident 2005
Robustness2_data$m <- robustness_2016sensis_df$FittedCurve
#Calcul du nombre de décès avec d = m * expo -> Exposition constante
Robustness2_data$d <- Robustness2_data$m * Robustness2_data$n

#Fit du modèle sur les nouvelles données 
SSE_robustness2<- morthump(data=Robustness2_data, model='sse')
#Récupération des coefficients alpha
SSE_coeffcients_robustness2 <-as.data.frame(as.matrix(coef(SSE_robustness2)))
#Récupération des taux fittés
SSE_deathrates_robustness2 <- as.data.frame(SSE_robustness2$mhat$mhat3[,1]+SSE_robustness2$mhat$mhat1[,1]+SSE_robustness2$mhat$mhat2[,1])
row.names(SSE_deathrates_robustness2) <- 0:85

#Plot des composantes et de la courbé fittée
hump_r2 <- melt(exp(SSE_robustness2$XX$X3%*%(SSE_coeffcients_robustness2[23:42,])),varnames = c('Age','Year'),value.name = 'm')
hump_r2$Comp <- 'Hump'
senescent_r2 <- melt(exp(SSE_robustness2$XX$X2%*%(SSE_coeffcients_robustness2[3:22,])),varnames = c('Age','Year'),value.name = 'm')
senescent_r2$Comp <- 'Senescent'
infant_r2<-  melt(exp(SSE_robustness2$XX$X1%*%(SSE_coeffcients_robustness2[1:2,])), varnames = c('Age','Year'),value.name = 'm')
infant_r2$Comp <- 'Infant'
#PLot du comportement des Qx Fittés
MU_r2 <- hump_r2
MU_r2$m <- hump_r2$m+infant_r2$m+senescent_r2$m
MU_r2$Comp <-'FittedCurve'

#Plot de toutes les composantes + des qx fittés
ALL_plot_robust2 <- rbind(hump_r2, infant_r2, senescent_r2, MU_r2)
ALL_plot_robust2$Year <-'2016_newfit'
ALL_plot_robust2$Age <- 0:85

ALL_plot_robust22 <- melt(robustness_2016sensis_df[,-6], id.vars = c('Age'), value.name = 'm', varnames = c('Infant','Senescent', 'Hump', 'FittedCurve'), variable.name = 'Comp')
ALL_plot_robust22$Year <- '2016_before'

ALL_plot_robust122 <- rbind(ALL_plot_robust2, ALL_plot_robust22)

plot_robust2 <- ggplot(subset(ALL_plot_robust122)) + 
  geom_line( aes(x=Age, y=m, group= interaction(Year,Comp), color=Comp,linetype=Year))+
  scale_y_continuous(trans='log10',limits = c(1e-5,1),name = 'μ')  +ggtitle('SSE Composantes')

#EV du nouveau fit 
LE_period_Model(Age = 0, Period = 'm', QxModel = MU_r2 )


#Graphique du fit 2004 avec la composante accident 2005
png('robust2MAR.png', width= 3500, height=2000, res=400)
plot_robust2
dev.off()
