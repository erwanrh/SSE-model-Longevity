######## SCRIPT MODEL SSE  #######-
# Author: Erwan Rahis (erwan.rahis@axa.com), GRM Life, Longevity Team
# Version 1.0
# Last update: 12/05/2020

# Script to plot life expectancies and raw data for SSE
# Camarda (2016)




####################################################################################-
#-------------------------- Raw Data Plot ----------------------------
####################################################################################-


# Plot des données brutes ------------------------------
All_data_male <- data.frame()
All_data_female<- data.frame()
for (year in 2000:2016){
  Data_males_temp <- HMD2MH(country=country_code,year=year, sex='males',path='Data',xtra=TRUE)
  Data_males_temp$year <- year
  All_data_male <- rbind(All_data_male, Data_males_temp)
  
  Data_females_temp <- HMD2MH(country=country_code,year=year, sex='females',path='Data',xtra=TRUE)
  Data_females_temp$year <- year
  All_data_female <- rbind(All_data_female, Data_females_temp)
}

All_data_male$sex <- 'M'
All_data_female$sex <- 'F'

All_data_both <- rbind(All_data_female, All_data_male)

dplot <- ggplot(All_data_both) + geom_line(aes(x= x, y=d, group = year, color = year)) + ggtitle('Décès') + facet_wrap(~ sex ) +
  xlab('')
nplot <- ggplot(All_data_both) + geom_line(aes(x= x, y=n, group = year, color = year)) + ggtitle('Expositions') + facet_wrap(~sex)+
  xlab('age')
mplot <- ggplot(All_data_both) + geom_line(aes(x= x, y=m, group = year, color = year), na.rm = T) + ggtitle('Taux de mortalité') + facet_wrap(~sex)+
  xlab('age')+ scale_y_continuous(trans = 'log10') +
  annotate("rect", xmin=13, xmax=30, ymin=5e-5 , ymax=4e-03, alpha=0.5,  fill="#F7DC6F")+
  annotate("rect", xmin=-1, xmax=10, ymin=2e-05, ymax= 5e-02, alpha=0.5, fill='#BB8FCE')


mplot
grid.arrange(dplot, nplot)

ratedhmd <- dcast(subset(All_data_both, sex=='M'), formula = x ~ year, value.var = 'm')
ratedhmd[ratedhmd == 0] <- NA
All_data_both[is.na(All_data_both$m) |  All_data_both$m== 0, 'm'] <- 1



####################################################################################-
#-------------------------- Life Expectancy Plots  ----------------------------
####################################################################################-

# Tableau Latex d'espérances de vie ------------------------------
LE_table_latex2 <- cbind(rbind(LE_period_Model(Age = 0, Period = 2000,SSE_deathrates_male_df),
                               LE_period_Model(Age = 0, Period = 2016, SSE_deathrates_male_df)), 
                         rbind(LE_period_Model(Age = 0, Period = 2000, SSE_deathrates_female_df), 
                               LE_period_Model(Age = 0, Period = 2016,SSE_deathrates_female_df)))

colnames(LE_table_latex2) <- c('Male SSE', 'Female SSE')
row.names(LE_table_latex2) <- c('2000', '2016')



library(xtable)
xtable(LE_table_latex2, type = "latex", file = "filename2.tex")


# Données d'espérance de vie et améliorations ------------------------------

ex_m <- data.frame()
ex_f <- data.frame()
for (year_ in 2000:2016){
  ex_m[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_male_df) 
  ex_f[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_female_df)
  
}

ex_m$var <- c(NA, diff(ex_m$V1))*12
ex_m$sex <- 'Male'
ex_m$year <- 2000:2016


ex_f$var <- c(NA, diff(ex_f$V1))*12
ex_f$sex <- 'Female'
ex_f$year <- 2000:2016


ex_tot <- rbind(ex_m, ex_f)


#Plot de l'ex totale avec gap HOMMES FEMMES ---------------
ex_plot <- ggplot(ex_tot) + geom_line(aes(x= year, y=V1, group = sex, color =sex)) + ylab('Ex') +
  annotate("segment", y = ex_m['2000', 'V1'], yend = ex_f['2000', 'V1'], 
           x = 2000, xend = 2000) +
  annotate("segment", y = ex_m['2016', 'V1'], yend = ex_f['2016', 'V1'], 
           x = 2016, xend = 2016) 

png('ex_deltaex_FR.png', width =5000, height = 3000, res =650 )
ex_plot
dev.off()


#Tableau de Moyenne des améliorations en ex (latex export) --------------
ex_means <- data.frame(mean(ex_m$var, na.rm = T), mean(ex_f$var, na.rm = T), row.names = 'Mean LE Improvement' )
colnames(ex_means) <- c('Males','Females')

xtable(ex_means, type = "latex", file = "filename2.tex")


#Plot des améliorations annuelles de l'ex HOMMES ET FEMMES  ---------------
ex_imp_plot <- ggplot(ex_tot, aes(x= year, y = var, group = sex, color = sex)) +
  geom_line() + ylab('LE Improvement (months)') +
  annotation_custom(tableGrob(round(ex_means, 2), cols = c('Male', 'Female'), theme=ttheme_default(base_size = 8)),
                    xmax=2018, ymin=4, ymax=4)+
  ggtitle('LE Improvements')

png('ex_imp.png', width =4000, height = 2500, res =550 )
ex_imp_plot
dev.off()


#Moyennes  sur 5 ans des améliorations en ex --------------- 
ex_imp_5mean_f <- aggregate(ex_f[-1,], list(rep(1:(nrow(ex_f)%/%5+1),each=5,len=nrow(ex_f)-1)),mean, na.action = na.omit)[-1][,c('var','year')]
ex_imp_5mean_f$sex <- 'Female'

ex_imp_5mean_m  <- aggregate(ex_m[-1,], list(rep(1:(nrow(ex_m)%/%5+1),each=5,len=nrow(ex_m)-1)),mean, na.action = na.omit)[-1][,c('var','year')]
ex_imp_5mean_m$sex <- 'Male'

ex_imp_5mean_tot <- rbind(ex_imp_5mean_f, ex_imp_5mean_m)


# PLot des moyennes sur 5 ans --------------
ex_imp_plot_5y <- ggplot() +
  geom_line(data = ex_imp_5mean_tot, aes(x= year, y = var, group = sex, color = sex)) +
  geom_point(data= ex_tot, aes(x= year, y = var, group = sex, color = sex))+
  ylab('LE Improvement (months)') + 
  ggtitle('Le Improvements and 5Y AVG')


png('ex_imp.png', width =4000, height = 4000, res =550 )
grid.arrange(ex_imp_plot, ex_imp_plot_5y)
dev.off()


#Moyennes mobiles  sur 5 ans des améliorations en ex --------------
ex_mean_rolling5 <- data.frame()
i<-1
for (year in 2000:2017){
  ex_mean_rolling5[i,'var'] <- mean(ex_f[ex_f$year <= year & ex_f$year > year-5, 'var'])
  ex_mean_rolling5[i+1,'var'] <- mean(ex_m[ex_m$year <= year & ex_m$year > year-5, 'var'])
  ex_mean_rolling5[i,'year'] <- year
  ex_mean_rolling5[i+1,'year'] <- year
  ex_mean_rolling5[i,'sex'] <- 'Female'
  ex_mean_rolling5[i+1,'sex'] <- 'Male'
  i <- i+2
}


# PLot des moyennes mobiles sur 5 ans--------------
ex_imp_plot_5yrolling <- ggplot() +
  geom_line(data = ex_mean_rolling5, aes(x= year, y = var, group = sex,linetype= sex), color = '#5DADE2') +
  geom_point(data= ex_tot, aes(x= year, y = var, group = sex, shape= sex), alpha = 0.5, )+
  scale_shape_manual('Improvements', values = c(16,16)) +
  scale_linetype_manual('Rolling 5Y Average', values = c(1,1)) +
  scale_x_continuous(breaks = seq(2000,2016,10), minor_breaks =seq(2000,2016,5))+
  ylab('LE Improvement (months)') + 
  ggtitle('LE Improvements and 5Y rolling AVG') +
  facet_wrap(~ sex, nrow = 2)


png('ex_rollingimp.png', width =5000, height = 4000, res =550 )
ex_imp_plot_5yrolling
dev.off()
