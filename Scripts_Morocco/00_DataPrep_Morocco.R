####################################################################################-
#-------------------------------------DATA PREP MOROCCO ------------------------------------
####################################################################################-

# Data qx ------------
#importing  interpolated data for males and females + importing Raw Data from WHO (5years age groups + replacing Infant <1yo mortality with raw data)
Data_qx_MAR_F <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_mx_polyfit_F.csv', row.names = 1)
Data_qx_MAR_F$sex <- 'F'
Data_qx_MAR_M <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_mx_polyfit_M.csv', row.names = 1)
Data_qx_MAR_M$sex <- 'M'
All_raw_data <- read.csv2('Data/MAR/alldata_MAR.csv',fileEncoding="UTF-8-BOM" , dec ='.')
All_raw_data$type <- 'Raw WHO Data'
Data_qx_MAR_F[Data_qx_MAR_F$Age == 0,'qx_male'] <- All_raw_data[((All_raw_data$age == 0) & (All_raw_data$sex=='F')),'mx']
Data_qx_MAR_M[Data_qx_MAR_M$Age == 0,'qx_male'] <- All_raw_data[((All_raw_data$age == 0) & (All_raw_data$sex=='M')),'mx']
Data_qx_MAR_All <- rbind(Data_qx_MAR_F, Data_qx_MAR_M)
Data_qx_MAR_All$type <- 'Interpolated Data'



# Plot Qx Raw
MAR_qx_plot <- ggplot(Data_qx_MAR_All) + 
  geom_line(aes(x= Age, y = qx_male, group=year, color=year, linetype = type)) + 
  geom_point(data = All_raw_data, aes(x=age, y= mx, group = year, color=year, shape = type))+
  scale_y_continuous(trans = 'log10') + 
  facet_wrap(~sex) + 
  ylab('qx') + ggtitle('mx Interpolated for Morocco - Polynomial Fit') +
  scale_x_continuous(breaks=seq(0,110, 10)) +
  scale_linetype_manual('', values = c(1)) +
  scale_shape_manual('', values = c(4)) +
  guides(colour = guide_colourbar(order = 3),shape = guide_legend(order = 2), linetype = guide_legend(order = 1))

png('MAR_qx_interpolated_plot.png', height = 1000, width = 3000, res = 200)
MAR_qx_plot
dev.off()




# Data dx --------
Data_dx_MAR_M <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_deaths_splines_M.csv', row.names = 1)
Data_dx_MAR_M$sex <- 'M'
Data_dx_MAR_F <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_deaths_splines_F.csv', row.names = 1)
Data_dx_MAR_F$sex <- 'F'
Data_dx_MAR_F[Data_dx_MAR_F$Age<= 10,'qx_male'] <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_deaths_polyfit_F.csv', row.names = 1)[Data_dx_MAR_F$Age<= 10,'qx_male']
Data_dx_MAR_M[Data_dx_MAR_M$Age <= 10,'qx_male'] <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_deaths_polyfit_M.csv', row.names = 1)[Data_dx_MAR_M$Age<= 10,'qx_male']
Data_dx_MAR_All <- rbind(Data_dx_MAR_F, Data_dx_MAR_M)
Data_dx_MAR_All$type <- 'Interpolated Data'


# Plot dx
MAR_dx_plot <- ggplot(Data_dx_MAR_All) + 
  geom_line(aes(x= Age, y = qx_male, group=year, color=year, linetype = type)) + 
  geom_point(data = All_raw_data, aes(x=age, y= deaths, group = year, color=year, shape = type))+ facet_wrap(~sex) + 
  ylab('qx') + ggtitle('Death count Interpolated for Morocco - Cubic Splines Fit') + #scale_colour_gradientn(colors = rainbow(5)) +
  scale_x_continuous(breaks=seq(0,110, 10)) +
  scale_linetype_manual('', values = c(1)) +
  scale_shape_manual('', values = c(4)) +
  guides(colour = guide_colourbar(order = 3),shape = guide_legend(order = 2), linetype = guide_legend(order = 1))

png('MAR_dx_interpolated_plot.png', height = 1000, width = 3000, res = 200)
MAR_dx_plot
dev.off()



# Data nx -------
Data_nx_MAR_M <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_expo_splines_M.csv', row.names = 1)
Data_nx_MAR_M$sex <- 'M'
Data_nx_MAR_F <- read.csv('Data/MAR/WHO_DATA/Morocco_plot_expo_splines_F.csv', row.names = 1)
Data_nx_MAR_F$sex <- 'F'
Data_nx_MAR_All <- rbind(Data_nx_MAR_F, Data_nx_MAR_M)
Data_nx_MAR_All$type <- 'Interpolated Data'

# Plot nx
MAR_nx_plot <- ggplot(Data_nx_MAR_All) + 
  geom_line(aes(x= Age, y = qx_male, group=year, color=year, linetype = type)) + 
  geom_point(data = All_raw_data, aes(x=age, y= expo, group = year, color=year, shape = type))+ facet_wrap(~sex) + 
  ylab('qx') + ggtitle('Expositions Interpolated for Morocco - Cubic Splines Fit') +  #scale_colour_gradientn(colors = rainbow(5)) +
  scale_x_continuous(breaks=seq(0,110, 10)) +
  scale_linetype_manual('', values = c(1)) +
  scale_shape_manual('', values = c(4)) +
  guides(colour = guide_colourbar(order = 3),shape = guide_legend(order = 2), linetype = guide_legend(order = 1))

png('MAR_nx_interpolated_plot.png', height = 1000, width = 3000, res = 200)
MAR_nx_plot
dev.off()

#Final Plot

png('MAR_RAWDATA_interpolated_plot.png', height = 2000, width = 3000, res = 300)
grid.arrange(MAR_dx_plot, MAR_nx_plot)
dev.off()


