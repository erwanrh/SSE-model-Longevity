
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 15:05:03 2020

@author: R-Erwan
"""



####Interpolation TEST
import pandas as pd
import numpy as np
import matplotlib.pyplot as mp
from scipy.interpolate import interp1d


#Data
path = 'C:/Users/r-erwan/Downloads/SSE-model-Longevity/'
Morocco_Deaths_7_all = pd.read_csv(path + 'Data/MAR/alldata_MAR.csv',  encoding = 'utf-8-sig'  ,sep=';')

#Variable to interpolate
interpolation_var = 'mx'
#Choice of sex
sex = 'F'
Morocco_Deaths_7 = Morocco_Deaths_7_all[Morocco_Deaths_7_all['sex']==sex]

#PLOT FUNCTION
def plotPolyFit(PolyFit):
    t= np.linspace(PolyFit['ad'],PolyFit['af'],PolyFit['af']-PolyFit['ad'])
    mp.plot(t,np.exp(PolyFit['PF'](t)),'.',PolyFit['X'],PolyFit['Y'],'o')
    mp.title('Morrocan Mortality ages:'+str(PolyFit['ad'])+'-'+str(PolyFit['af'])+' year:'+str(PolyFit['YR']))


#FONCTION POLYFIT
def Polyfit_(ad,af,YR,deg):
    selection= Morocco_Deaths_7[(Morocco_Deaths_7['year']==YR) & (Morocco_Deaths_7['age']>=ad) &   
                                (Morocco_Deaths_7['age']<=af)].sort_values(by='age')
    Y=selection[interpolation_var]
    X=selection['age']
    PF=np.polynomial.polynomial.Polynomial.fit(x=X,y=np.log(Y), deg=deg)
    res=np.exp(PF(np.linspace(ad,af,af
                              
                              -ad+1)))
    return  {'PF':PF,'X':X,'Y':Y,'ad':ad,'af':af,'YR':YR,'res':res}
    

#FONCTIONS FIT
def Extrapolation_Fit(A1,A2,A3,YR, save, plot):
    #AGES JEUNES
    PolyFit1=Polyfit_(0,A1,YR,3)
    
    #plotPolyFit(PolyFit1)
    
    PolyFit2=Polyfit_(A1+1,A2,YR,2)
    #plotPolyFit(PolyFit2)
    
    PolyFit3=Polyfit_(A2+1,A3,YR,2)
    #plotPolyFit(PolyFit3)
    #PolyFit4=Polyfit_(A3,110,YR,2)
    
    Interp_results = {}
    Interp_results[0]= PolyFit1['res']
    Interp_results[A1]= PolyFit2['res']
    Interp_results[A2]= PolyFit3['res']
    #Interp_results[A3]= PolyFit4['res']
    
    #ALL AGES OBSERVED DATA
    interpolation_data = Morocco_Deaths_7[(Morocco_Deaths_7['year']==YR)][interpolation_var]
    
    #ALL AGES FY INTERPOLATED DATA
    Interp_FY = np.concatenate([Interp_results[0],Interp_results[A1],Interp_results[A2]]) #, Interp_results[A3]] )
    
    #TEST : Spline
    splines3d = interp1d(x=Morocco_Deaths_7[(Morocco_Deaths_7['year']==YR)]['age'],y=interpolation_data,kind='cubic')
    
    ##ALL PLOTS SPLINES + POLYFIT
    if plot == True:      
        mp.clf()
        mp.plot(np.linspace(0,110,111),splines3d(np.linspace(0,110,111)),color='#3498DB')
        mp.plot(np.linspace(0,85,86),Interp_FY,color='#AF7AC5')
        mp.plot( Morocco_Deaths_7[(Morocco_Deaths_7['year']==YR)]['age'],interpolation_data,'.', alpha=0.7,color='#00A658')
        mp.legend(['SplinesCubiques','PolyFit','ObsData'],loc='upper left')
        mp.grid(color='grey',linestyle='-',linewidth=0.1)
        if interpolation_var == 'qx' or  interpolation_var == 'mx' :
            mp.yscale('log')
        mp.title('Morrocan Mortality Segments:'+str(A1)+'-'+str(A2)+' year:'+str(YR))
        if save==True:
            mp.savefig(path + 'Data/MAR/Interpolate_test_'+ interpolation_var + '_'+str(A1)+'_'+str(A2)+'_'+str(A3)+'_'+str(YR)+'.png',format='png',dpi=500)
    
    return Interp_FY
   
    

 
 
def Splines_Fit(YR):
    interpolation_data = Morocco_Deaths_7[(Morocco_Deaths_7['year']==YR)][interpolation_var]
    splines3d = interp1d(x=Morocco_Deaths_7[(Morocco_Deaths_7['year']==YR )]['age'],y=interpolation_data,kind='cubic')
    Data=splines3d(np.linspace(0,110,111))
    return Data    




#DECOUPAGE DES AGES:
DataTable_Interpolated = pd.DataFrame()
DataTable_Splines = pd.DataFrame()


for year in range(2000, 2005):
    DataTable_Interpolated[year]=Extrapolation_Fit(20,40,85,year, False, False)
    DataTable_Splines[year]=Splines_Fit(year)


DataTable_Interpolated[2005]=Extrapolation_Fit(20,50,85,2005, False, True)

for year in range(2006, 2017):
    DataTable_Interpolated[year]=Extrapolation_Fit(20,58,85,year, False,False)
    DataTable_Splines[year]=Splines_Fit(year)


Extrapolation_Fit(20,50,85,2000, True, True)


#Ploynomial Fit
DataTable_Interpolated.to_csv(path + 'Data/MAR/WHO_DATA/Morocco_' + interpolation_var + '_polyfit_' + sex + '.csv')
Data_melted= pd.melt(DataTable_Interpolated.reset_index(),id_vars='index', value_name ='qx_male',var_name='year').rename(columns={'index':'Age'})
Data_melted.to_csv(path + 'Data/MAR/WHO_DATA/Morocco_plot_'+ interpolation_var + '_polyfit_' + sex + '.csv')

#SPLINES
DataTable_Splines.to_csv(path + 'Data/MAR/WHO_DATA/Morocco_' + interpolation_var + '_splines_' + sex + '.csv')
Data_melted_splines= pd.melt(DataTable_Splines.reset_index(),id_vars='index', value_name ='qx_male', var_name='year').rename(columns={'index':'Age'})
Data_melted_splines.to_csv(path + 'Data/MAR/WHO_DATA/Morocco_plot_'+ interpolation_var + '_splines_' + sex + '.csv')
