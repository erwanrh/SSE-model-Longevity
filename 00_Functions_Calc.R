#Compute Life Expectancies from Qx
QX_model_period <- function(Age,Period,QxModel){  #Calcul des Qx en cohorte pour un age donn? ------------------------------------------------------------------------
  
  Qx_Period<- QxModel[as.numeric(row.names(QxModel))>=Age, colnames(QxModel) == Period]
  Qx_Period
}


tpx_period <- function(Age,Period, QxModel){ #Calc tpx cohort from qx table ------------------------------------------------------------------------
  Qx_Period <- QX_model_period(Age,Period,QxModel)
  tpx<-data.frame()
  tpx[1,1]<-1-Qx_Period[1]
  
  for (i in 2:length(Qx_Period)){
    tpx[1,i]<-tpx[1,i-1]*(1-Qx_Period[i])
  }
  
  tpx_Cohort<- tpx
}



LE_period_Model <- function(Age,Period, QxModel){ # LE cohort to compare to compute indicators ------------------------------------------------------------------------
  sum(tpx_period(Age,Period, QxModel),na.rm = TRUE ) +0.5
}





#Function to compute entropy on period 
compute_entropy_period <- function(Qx.Table, start.age)
{
  entropy_df <- data.frame()
  
  for (age in start.age){
    lx <- compute_lx_period(Qx.Table, age)
    ln_lx <- log(lx/lx[1,])
    
    h1 <- colSums(lx * ln_lx, na.rm = TRUE)
    h2 <- colSums(lx)
    
    
    entropy_df_temp <- data.frame(row.names = NULL, period = as.character(names(-h1/h2)), value=-h1/h2)
    entropy_df_temp$age <- age
    entropy_df <- rbind(entropy_df, entropy_df_temp)   
  }
  
 
  entropy_df$variable <- 'entropyPeriod'
  
  entropy_df
  
}


compute_lx_period <- function(Qx.Table, start.age) #Compute survival function over a period
{
  #Max age of Qx Table
  age_max <- as.numeric(tail(row.names(Qx.Table), n = 1)) 
  
  #Warning if the age given is superior to the max age in the qx table
  if (age_max  <= start.age){
    warning("Age given in parameters is superior or equal to maximum table age.")
    return(NaN)
  }
  
  px_table <- rbind(rep(1, ncol(Qx.Table)), (1 -Qx.Table[as.numeric(row.names(Qx.Table)) >= start.age,]) )
  
  apply(px_table, MARGIN = 2, FUN = cumprod)*100000
}

