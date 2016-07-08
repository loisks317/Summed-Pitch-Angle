# read_rbsp_output.py
#
#  generate the plots in a faster and more accesible fashion
import os
import glob 
import numpy as np
from matplotlib import pyplot as plt
import datetime
import matplotlib.dates as dates
import spacepy.pybats.kyoto as spk
from spacepy import pycdf
import itertools as itert
import math
from numpy import ma
import pandas as pd
# 
# SETTINGS
n1=0
n2=16
PAs=11 # 11 pitch angle bins on HOPE
dateStart='20130201'           # starting date
#dateEnd='20130301'
dateEnd='20150401'
#dateEnd='20150401'             # ending date
nLbins=11                       # for the number of L shells, spaced by .25
LbinsArr=np.linspace(1.5, 4,nLbins)
nmlt_bins=48                   # 30 minute time resolution in MLT
MLTbinsArr=np.linspace(0.25, 23.75, nmlt_bins)
satellite=['A', 'B']
name_species=[ 'FPDU', 'FHEDU', 'FODU']
spe=['P','He','O']
mass=[2, 4, 16]
#
# Data arrays
def weightedAvg(window,df,weight,DT):
   # window = time window
   # df = data frame to avg over
   # weight = array of weights
   # DT = starting date time
   cutoff=1440/window # number of windows in a day by minutes
   newDF={'MLT':np.zeros(cutoff),'L':np.zeros(cutoff), 'flux':np.zeros(cutoff), 'energy':np.zeros(cutoff), 'delta':np.zeros(cutoff) }
   weight=np.array(weight)
   for iTime in range( cutoff-1):
      sI=df.index.searchsorted(DT+iTime*datetime.timedelta(minutes=5))
      eI=df.index.searchsorted(DT+(iTime+1)*datetime.timedelta(minutes=5))
      # we have five minute window now get the weights
      totalC=1.0*np.nansum(weight[sI:eI])
      newDF['MLT'][iTime]=np.nansum(np.array(df['MLT'])[sI:eI]*weight[sI:eI])/totalC
      newDF['L'][iTime]=np.nansum(np.array(df['L'])[sI:eI]*weight[sI:eI])/totalC
      newDF['flux'][iTime]=np.nansum(np.array(df['flux'])[sI:eI]*weight[sI:eI])/totalC
      newDF['energy'][iTime]=np.nansum(np.array(df['energy'])[sI:eI]*weight[sI:eI])/totalC
      newDF['delta'][iTime]=np.nansum(np.array(df['delta'])[sI:eI]*weight[sI:eI])/totalC
   newDF['MLT'][-1]=df['MLT'][-1]
   newDF['L'][-1]=df['L'][-1]
   newDF['flux'][-1]=df['flux'][-1]
   newDF['energy'][-1]=df['energy'][-1]
   newDF['delta'][-1]=df['delta'][-1]
   timeArr=pd.date_range(DT,DT+datetime.timedelta(days=1), freq='5Min')[:-1]
   DF=pd.DataFrame(newDF, index=timeArr)
   # 5 minute weighted average
   return DF
# prep the start times
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
for iSpe in range(len(name_species)):
 dfs={} # for a and b data storage
 for iSat in range(len(satellite)):
  fplasmaDensity=[]
  fMLT=[]
  fL=[]
  fkp=[]
  DT=datetime.datetime.strptime(dateStart, '%Y%m%d') 
  while DT != endDt:
     date=datetime.datetime.strftime(DT, '%Y%m%d')
     DT=datetime.datetime.strptime(date, '%Y%m%d')


     # load in HOPE files
     os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_L3_'+satellite[iSat])
     f=glob.glob('*'+date+'*')
     try:
         pyf=pycdf.CDF(f[0])
         #
         # put into pandas arrays
         energy=np.swapaxes(pyf['HOPE_ENERGY_Ion'][...],1,0)
         epoch=pyf['Epoch_Ion'][...]
         MLT=pyf['MLT_Ion'][...]
         L=pyf['L_Ion'][...]
         EnergyDelta=np.swapaxes(pyf['ENERGY_Ion_DELTA'][...],1,0)[n1:n2]
         FLUX=pyf[name_species[iSpe]][...][:,5,:] # get the PA = 90
         FLUXl=np.swapaxes(FLUX,1,0)[n1:n2]
         dataFrames={}
         kyoto=spk.fetch('kp', DT, DT)
         day=datetime.datetime.strftime(DT,'%d') # need to get day of month
         kp_arr=kyoto['kp'][(int(day)-1)*8: int(day)*8]
         kpC=0
         kpM=0
         nKP=np.zeros(288)
         for iKP in range(288):
              nKP[iKP]=kp_arr[kpM]
              kpC+=1
              if kpC == 36:
                 kpC=0
                 kpM+=1
              else:
                 continue

         rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
         rng=rt[::1]
         kyotodf2=pd.DataFrame({'kp':nKP}, index=rng)
         #kyotodf2=kyotodf['kp'].resample('5min',how='median').reindex(index=rng)
         COUNTS=np.swapaxes(pyf['Counts_'+spe[iSpe]+'_Omni'][...],1,0)[n1:n2]
         #
         # each dataFrame is set up to be by the energy channel number (so 0 = 0.9 eV)
         for iN in range(n1,n2):
            index=np.where(energy[iN] > 20)[0]
            energy[iN][index]=np.nan
            FLUXl[iN][index]=np.nan # bad energies
            COUNTS[iN][index]=np.nan # bad energies
            offFlux=np.where((FLUXl[iN] <= 0) & (FLUXl > 1e12))[0]
            FLUXl[iN][offFlux]=np.nan
            COUNTS[iN][offFlux]=np.nan
            dataFrames[str(iN)]=pd.DataFrame({'MLT':MLT, 'L':L, 'flux':FLUXl[iN], 'energy':energy[iN], 'delta':EnergyDelta[iN]}, index=epoch)
            # now resample the data frame
            #
            # NEED TO WRITE MY OWN CODE TO RESAMPLE FOR 5 MINUTES BUT FACTOR IN WEIGHTED AVERAGE

            dataFrames[str(iN)]=weightedAvg(5, dataFrames[str(iN)], COUNTS[iN], DT)
            dataFrames[str(iN)]['kp']=kyotodf2 # add resampled kp onto the resampled flux
            #
            # now sum for plasma density
            Energy_Ev=dataFrames[str(iN)]['energy']*(1.602*1e-19)
            mass_adj=mass[iSpe]*(1.67*1e-27) 
            dataFrames[str(iN)]['density']=4.0*np.pi*(1.0/np.sqrt(2.0*Energy_Ev/mass_adj))*np.array(dataFrames[str(iN)]['delta'])*(1.0e-5)*np.array(dataFrames[str(iN)]['flux'])
         #plasmaDensity=list(np.array(dataFrames['0']['density']))
         #MLT=list(np.array(dataFrames['0']['MLT']))
         #L=list(np.array(dataFrames['0']['L']))
         #kp=list(np.array(dataFrames['0']['kp']))
         for iN in range(n1,n2):
           fplasmaDensity+=list(np.array(dataFrames[str(iN)]['density']))
           fMLT+=list(np.array(dataFrames[str(iN)]['MLT']))
           fL+=list(np.array(dataFrames[str(iN)]['L']))
           fkp+=list(np.array(dataFrames[str(iN)]['kp']))
     except(OSError):
        print('bad date: ' + str(date))
     DT=DT+datetime.timedelta(days=1)
     #
     # now sort
  kpIndex=np.where(np.array(fkp)>=3)[0]
  try:
     fplasmaDensity=np.array(fplasmaDensity)
     fplasmaDensity[kpIndex]=np.nan
  except(TypeError):
     fplasmaDensity=fplasmaDensity
     print("error")
  sPden=[ [[] for i in range(nLbins)] for j in range(nmlt_bins)]
  for iL in range(nLbins-1):
       for iMLT in range(nmlt_bins-1):
         pL=np.array(fL); pMLT=np.array(fMLT)
         Lindex=np.where((pL>=LbinsArr[iL]-0.125) & (pL<LbinsArr[iL+1]-0.125))[0]
         MLTindex=np.where((pMLT>=MLTbinsArr[iMLT]-0.25) & (pMLT<MLTbinsArr[iMLT+1]-0.25))[0]
         # get the matches
         matches=list(set(Lindex) & set(MLTindex))
         # remove the nans
         try:
            sPden[iMLT][iL]=np.array(fplasmaDensity)[matches]#[~np.isnan(plasmaDensity[matches])]
            sPden[iMLT][iL]=sPden[iMLT][iL][~np.isnan(sPden[iMLT][iL])]
         except(TypeError):
            sPden[iMLT][iL]=[]
  dfs[satellite[iSat]]=sPden
     # save this file
  os.chdir('/Users/loisks/Documents/Functions/')
  import pickling
  os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
  pickling.hdf5_data_save(sPden, 'PA=90_plasma_density_'+str(n1)+'_'+str(n2)+'_sat='+satellite[iSat]+'_species='+spe[iSpe], 'PlasmaDensity', nmlt_bins, nLbins)
  # combined files
 os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
 cDen=[ [[] for i in range(nLbins)] for j in range(nmlt_bins)]
 for iL in range(nLbins-1):
   for iMLT in range(nmlt_bins-1):
      cDen[iMLT][iL]+=list(dfs['A'][iMLT][iL])+list(dfs['B'][iMLT][iL])
 pickling.hdf5_data_save(cDen, 'PA=90_combined_plasma_density_'+str(n1)+'_'+str(n2)+'_species='+spe[iSpe], 'PlasmaDensity', nmlt_bins, nLbins)      
         
         

 












## CHECK OVER THIS CODE IN THE MORNING
