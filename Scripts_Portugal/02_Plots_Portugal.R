######## SCRIPT MODEL SSE Plot Functions #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to fit a Sums Of Smooth Exponential model according to the paper "Sum Of Smooth Exponentials to decompose complex series of counts" 
# Camarda (2016)




####################################################################################-
#-------------------------- Plot funtions MOROCCO  ----------------------------
####################################################################################-



# 3.i Plot des courbes de mortalité  ----------------------------
DeathRates_plot_males <- melt(SSE_deathrates_male_df,varnames = c('Age','Year'),value.name = 'coeff',)
DeathRates_plot_females <- melt(SSE_deathrates_female_df,varnames = c('Age','Year'),value.name = 'coeff',)
DeathRates_plot_males$gender <- 'M'
DeathRates_plot_females$gender <- 'F'
DeathRates_plot <- rbind(DeathRates_plot_females,DeathRates_plot_males)
DeathRates_plot$Age <- as.numeric(DeathRates_plot$Age)

DeathRatesplot <- ggplot(DeathRates_plot)+geom_line(aes(x=Age, y=coeff, group = Year,color=Year))+ ggtitle(paste0('SSE Mortality Rates - ', country_code)) +facet_wrap(gender~.,scales = 'free')+
  scale_y_continuous(trans='log2') + scale_x_continuous(breaks = seq(0,110,5), labels= seq(0,110,5) )

ggsave(paste0('drplot_',country_code,'.pdf'),DeathRatesplot,width = 20, height = 12)
print('Qx Plot Saved', quote= FALSE)

# 3.ii Plot des coefficients  ----------------------------
COMP1 <- paste0('X',1:2)
COMP2 <- paste0('X',3:17)
COMP3 <- paste0('X',18:42)

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
print('Coefficients Plot Saved', quote= FALSE)



# 3.iv Graph des élements des B-splines  ----------------------------
splines_plot <- SSE_splines_male
splines_plot [splines_plot$splinenb %in% 1:2, 'comp']<-'c1'
splines_plot[splines_plot$splinenb %in% 3:17, 'comp']<-'c2'
splines_plot[splines_plot$splinenb %in% 18:42, 'comp']<-'c3'

splines_ggplot <- ggplot(splines_plot,aes(x=age, y=splinevalue, color=year, group=year))+ geom_line() + facet_wrap(.~splinenb+ comp) + ggtitle(paste0('Splines - ', country_code))
ggsave(paste0('spline_comp_',country_code,'.pdf'),plot = splines_ggplot,width = 15, height = 15)


# 3.v Graph de tous les coefficients pour comparaison avec les splines  ----------------------------
coefs_plots_52gridM <- ggplot(Coeff_plot_males,aes(x=Year, y=coeff, color=year, group=year))+ geom_line() + facet_wrap(.~Coeff.num) + ggtitle (paste0('Coefficients for splines Males - ', country_code))
ggsave(paste0('coeffs_52grid_M',country_code,'.pdf'),plot = coefs_plots_52gridM,width = 15, height = 15) 

coefs_plots_52gridF <- ggplot(Coeff_plot_females,aes(x=Year, y=coeff, color=year, group=year))+ geom_line() + facet_wrap(.~Coeff.num) + ggtitle (paste0('Coefficients for splines Females - ', country_code))
ggsave(paste0('coeffs_52grid_F',country_code,'.pdf'),plot = coefs_plots_52gridF,width = 15, height = 15) 

print('Grid plots for B-Splines and Coefficients saved')