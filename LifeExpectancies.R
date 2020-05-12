#Compute Life Expectancies from Qx


QX_model_period <- function(Age,Period,QxModel){  #Calcul des Qx en cohorte pour un age donn? ------------------------------------------------------------------------
  
  Qx_Period<- QxModel[as.numeric(row.names(QxModel))>=Age,colnames(QxModel)==Period]
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
  sum(tpx_period(Age,Period, QxModel),na.rm = TRUE )
}





