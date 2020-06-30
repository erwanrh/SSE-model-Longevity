# Plot des données ------------------------------
All_data_male <- data.frame()
All_data_female<- data.frame()
for (year in 1960:2017){
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
# Tableau d'espérances de vie ------------------------------
LE_table_latex1 <- cbind(rbind(LE_period_Model(Age = 0, Period = 1960, dcast(All_data_male,value.var = 'm', x ~ year)),
                               LE_period_Model(Age = 0, Period = 2017, dcast(All_data_male,value.var = 'm', x ~ year))), 
                         rbind(LE_period_Model(Age = 0, Period = 1960, dcast(All_data_female,value.var = 'm', x ~ year)), 
                               LE_period_Model(Age = 0, Period = 2017, dcast(All_data_female,value.var = 'm', x ~ year))))


colnames(LE_table_latex1) <- c('Male HMD', 'Female')
row.names(LE_table_latex1) <- c('1960', '2017')



LE_table_latex2 <- cbind(rbind(LE_period_Model(Age = 0, Period = 1960,SSE_deathrates_male_df),
                              LE_period_Model(Age = 0, Period = 2017, SSE_deathrates_male_df)), 
                        rbind(LE_period_Model(Age = 0, Period = 1960, SSE_deathrates_female_df), 
                              LE_period_Model(Age = 0, Period = 2017,SSE_deathrates_female_df)))
                        
colnames(LE_table_latex2) <- c('Male SSE', 'Female SSE')
row.names(LE_table_latex2) <- c('1960', '2017')


LE_table_latex <- cbind(LE_table_latex1,LE_table_latex2)

library(xtable)
xtable(LE_table_latex, type = "latex", file = "filename2.tex")


# Gap d'espérance de vie------------------------------

gap_m <- data.frame()
gap_f <- data.frame()
for (year_ in 1960:2017){
  gap_m[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_male_df) - LE_period_Model(Age = 0, Period = year_, dcast(All_data_male,value.var = 'm', x ~ year))
  gap_f[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_female_df) - LE_period_Model(Age = 0, Period = year_, dcast(All_data_female,value.var = 'm', x ~ year))
  
}

gap_m$sex <- 'Male'
gap_m$year <- 1960:2017
gap_f$year <- 1960:2017
gap_f$sex <- 'Female'
gap <- rbind(gap_m, gap_f)

gap_plot <- ggplot(gap) + geom_line(aes(x= year, y=V1, group = sex, color =sex)) + ylab('Delta Ex') +
  annotate("text", x = 2000, y = 0.01, label = paste0("Mean F = ", round(mean(gap_f$V1),3) )) +
  annotate("text", x = 2000, y = 0.009, label = paste0("Mean M = ", round(mean(gap_m$V1),3) ))


# Plot d'espérance de vie------------------------------

ex_m <- data.frame()
ex_f <- data.frame()
for (year_ in 1960:2017){
  ex_m[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_male_df) 
  ex_f[as.character(year_),1] <- LE_period_Model(0, year_, SSE_deathrates_female_df)
  
}

ex_m$sex <- 'Male'
ex_f$year <- 1960:2017
ex_m$year <- 1960:2017
ex_f$sex <- 'Female'
ex_tot <- rbind(ex_m, ex_f)

ex_plot <- ggplot(ex_tot) + geom_line(aes(x= year, y=V1, group = sex, color =sex)) + ylab('Ex') +
  theme(legend.position = "none")+
  annotate("segment", y = ex_m['1960', 'V1'], yend = ex_f['1960', 'V1'], 
           x = 1960, xend = 1960) +
  annotate("segment", y = ex_m['2017', 'V1'], yend = ex_f['2017', 'V1'], 
           x = 2017, xend = 2017) 

grid.arrange(ex_plot, gap_plot, nrow=1)
