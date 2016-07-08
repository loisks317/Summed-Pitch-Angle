# PACorrectedSort.py
#
# Sort by MLT, L-SHELL, and Pitch angle
# this will, unfortunately, take a while :(
# correct for s/c potential
# 1 minute sum
#
# LKS November 2015. 
#
# imports 
import numpy as np
from spacepy import pycdf
import glob
import os
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
import pickle 
import h5py
import pandas as pd

os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/PA_check')
#
# parameters
date1='20130201'
date2='20150401'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
Lbins=11
mltbins=48
LbinsArr=np.linspace(1.5, 4.25,Lbins+1)
MLTbinsArr=np.linspace(0.25, 24.25,mltbins+1)
kpBar=1
# energy channels
n1=0
n2=16
PAs=11 # 11 pitch angles on HOPE
dir=['HOPE_PA_A', 'HOPE_PA_B']
name=['rbspa', 'rbspb']
satellite=['A','B']
species=[ 'P']#, 'He', 'O']
HopeEnergies=[.99,1.2,1.34,1.55,1.83,2.18,2.53,2.95,3.38,3.94,4.64,5.35,6.26,7.32,8.51,9.92]
name_species=['FPDU']#, 'FHEDU','FODU'] # gives the unidirectional flux with PA
PA_LABEL = [4.5] # this is over written later

# pickle function
#
def save_data_pickle(data, EN_LABEL, data_label) :
#
    import pickle
    import os
    with open(data_label + '_' +EN_LABEL, 'w') as f:
        pickle.dump(data, f)
    os.chdir('..')

                    
def weightedAvg(window,df,weight,DT, label):
      # window = time window
      # df = array to avg over
      # weight = array of weights
      # DT = starting date time
      # label = dataFrame label
      cutoff=1440/window # number of windows in a day by minutes
      newDF=np.zeros(cutoff)
      weight=np.array(weight)
      for iTime in range( cutoff-1):
         sI=df.index.searchsorted(DT+iTime*datetime.timedelta(minutes=1))
         eI=df.index.searchsorted(DT+(iTime+1)*datetime.timedelta(minutes=1))
         # we have five minute window now get the weights
         totalC=1.0*np.nansum(weight[sI:eI])
         newDF[iTime]=np.nansum(np.array(df[label])[sI:eI]*weight[sI:eI])/totalC
      newDF[-1]=df[label][-1]
      timeArr=pd.date_range(DT,DT+datetime.timedelta(days=1), freq='1Min')[:-1]
      #DF=pd.DataFrame(newDF, index=timeArr)
      # 1 minute weighted average
      return newDF
#
#
# Need to set up that it saves every month
# Save by pitch angle
# save by s/c corrected energy
#    
# energy channels
for ispe in range(len(name_species)):
  for idir in range(2):# correct this
   # start overarching loop through days
   PA_Flkp=[ [ [  [ [] for x in xrange(Lbins)] for y in xrange(mltbins)] for z in  range(n1,n2)] for w in range(PAs)] # Energy x PA x MLT x L
#   os.chdir('/Users/loisks/Desktop/PA_Check/CorrectedPAFluxes/')
#   check=glob.glob()
#   try:
#       check[0] # see if anything is there ,this will result in an index error if nothing is there
#       print species[ispe]+'_'+name[idir]+'_'+str(ien)
#   except(IndexError):
   dt1=dt0
   while dt1 < dt2:
        monthCur=dt1.month
        yrCur=dt1.year
        try:
                date1=datetime.datetime.strftime(dt1, '%Y%m%d')
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/'+dir[idir])
                hope_file=name[idir]+'_rel03_ect-hope-PA-L3_'+date1+'*.cdf'
                ghope=glob.glob(hope_file)
                pyf=pycdf.CDF(ghope[0])
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+satellite[idir])
                f=glob.glob('*'+date1+'*')
                pyfe=pycdf.CDF(f[0])

        # get the data quantities we need, separate by species 
                # now a 2 d array, chosen for correct PA
                
            # now resample and rebin
                epoch=pyf['Epoch_Ion'][...] 
                L=pd.DataFrame({'L':pyf['L_Ion'][...]}, index=epoch) # L Shell
                energy=pyf['HOPE_ENERGY_Ion'][...] # Hope Energy
                MLT=pd.DataFrame({'MLT':pyf['MLT_Ion'][...]}, index=epoch) # MLT
                   
                PA_LABEL=pyf['PITCH_ANGLE'][...] # right PA label
                
                #
                # resample by a minute
                rt = pd.period_range(date1,periods=1440, freq='T').to_timestamp()
                rng=rt[::1]
                dfFlux={}
                dfCounts={}
                dfEnergy={}
                wMLT={}; wL={}; wEnergy={}; wFlux={}
                counter1=0
                kyoto=spk.fetch('kp', dt1, dt1)
                day=datetime.datetime.strftime(dt1,'%d') # need to get day of month
                kp_arr=kyoto['kp'][(int(day)-1)*8: int(day)*8]
                kpC=0
                kpM=0
                nKP=np.zeros(1440)
                for iKP in range(1440):
                    nKP[iKP]=kp_arr[kpM]
                    kpC+=1
                    if kpC == 180:
                        kpC=0
                        kpM+=1
                    else:
                        continue
                #kpdf=np.array(kyotodf2['kp'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
                #print nKP
                highKp=np.array(np.where(np.array(nKP)>kpBar)[0])
                #print highKp
                #print 'length of high Kp:' + str(len(highKp))
                # goes by pitch angle and initial energy channel
               # for iPA in range(PAs):
               #     for iEn in range(n1,n2):
                wFlux=[[[] for i in range(PAs)] for j in range(n1,n2)]
                wEnergy=[[[] for i in range(PAs)] for j in range(n1,n2)]
                for iPA in range(PAs):
                   for iEn in range(n1,n2):
                        ftemp=np.swapaxes(np.array(pyf[name_species[ispe]][...]),2 ,0)[iEn][iPA]
                        ftemp[ftemp<0]=np.nan
                        ftemp[ftemp>1e12]=np.nan
                        ftemp[np.swapaxes(energy, 1,0)[iEn]>20]=np.nan
                        Etemp=np.swapaxes(energy,1,0)[iEn]
                        Etemp[Etemp>20]=np.nan
                        dfFlux[str(counter1)]=pd.DataFrame({'Flux':ftemp}, index=epoch)
                        
                        dfCounts[str(counter1)]=pd.DataFrame({'Counts':np.swapaxes(pyf['Counts_'+species[ispe]][...],2,0)[iEn][iPA]}, index=epoch)
                        dfEnergy[str(counter1)]=pd.DataFrame({'Energy':Etemp}, index=epoch)
                        wMLT[str(counter1)]=np.array(MLT.resample('1min', how='mean').reindex(rng, fill_value=np.nan))
                        wL[str(counter1)]=np.array(L.resample('1min', how='mean').reindex(rng, fill_value=np.nan))
                        wEnergy[iEn][iPA]=np.swapaxes(np.array(dfEnergy[str(counter1)].resample('1min', how='mean').reindex(rng, fill_value=np.nan)),1,0)[0]
                        wFlux[iEn][iPA]=np.swapaxes(np.array(dfFlux[str(counter1)].resample('1min', how='mean').reindex(rng, fill_value=np.nan)),1,0)[0]
                        #
                        # weighted averages
                        #wMLT[str(counter1)]=weightedAvg(1, MLT, dfCounts[str(counter1)]['Counts'], dt1, 'MLT')
                        #wL[str(counter1)]=weightedAvg(1,L,dfCounts[str(counter1)]['Counts'],dt1, 'L')
                        #wEnergy[str(counter1)]=weightedAvg(1,dfEnergy[str(counter1)],dfCounts[str(counter1)]['Counts'],dt1, 'Energy')
                        #wFlux[str(counter1)]=weightedAvg(1,dfFlux[str(counter1)],dfCounts[str(counter1)]['Counts'],dt1, 'Flux')
                        # nan the high Kps
                        wFlux[iEn][iPA][highKp]=np.nan
                        #
                        # now nan high kp times
                        counter1+=1
                
                #
                # now get the EFW data and sort
                potential=-1*pyfe['Vavg'][...]
                potential[potential<-50]=np.nan
                potential[potential>10]=np.nan
                pEpoch=pd.DatetimeIndex(pyfe['epoch'][...])
                df=pd.DataFrame(potential, index=pEpoch, columns=['potential'])
                Phi=np.array(df['potential'].resample('1min', how='mean').reindex(rng, fill_value=np.nan))
                Phi=np.absolute(Phi)
                # nan the times of significant s/c charging
                Phi[Phi>5]=np.nan
                counter2=0
                # now we need to rebin the fluxes based on s/c corrected energy
                FluxC=[ [ [] for i in range(n1,n2)] for j in range(PAs)]
                MLTC=[ [ [] for i in range(n1,n2)] for j in range(PAs)]
                LC=[ [ [] for i in range(n1,n2)] for j in range(PAs)]
                #
                # everything is on a 1440 scale now
                for iPA in range(PAs):
                    CorrectEn=[ [] for i in range(n1,n2)]
                    for iEn in range(n1,n2):
                      CorrectEn[iEn]=np.array(wEnergy[iEn][iPA])+np.array(Phi)
                    CorrectEn=np.swapaxes(CorrectEn, 1, 0)
                    for iEn in range(n1,n2):
                        for iTime in range(len(CorrectEn)):
                              try:
                                  EnArr=np.array(CorrectEn[iTime])[~np.isnan(CorrectEn[iTime])]
                                  below=np.max(EnArr[EnArr<=HopeEnergies[iEn]])
                                  BelowI=np.where(EnArr==below)[0]
                                  above=np.min(EnArr[EnArr>=HopeEnergies[iEn]])
                                  AboveI=np.where(EnArr==above)[0]
                                  f=(HopeEnergies[iEn]- below)/(1.0*((HopeEnergies[iEn]-below)+(above-HopeEnergies[iEn])))
                                  intpFlux=(wFlux[AboveI][iPA][iTime]**f)*(wFlux[BelowI][iPA][iTime]**(1-f))

                              except(ValueError): # when there is no flux b/c of s/c charging
                                  intpFlux=np.nan
                              #cEn=min(range(len(HopeEnergies)), key=lambda i: abs(HopeEnergies[i]-wEnergy[str(counter2)][iTime]))
                              
                              FluxC[iPA][iEn].append(intpFlux)
                              MLTC[iPA][iEn].append(wMLT[str(counter2)][iTime])
                              LC[iPA][iEn].append(wL[str(counter2)][iTime])
                            # this should give an index error occasionally
                        counter2+=1
                
                    
                for iPA in range(PAs):
                  for iEn in range(n1,n2):
                     for iMLT in range(mltbins):
                       for iL in range(Lbins):
                           Lindex=np.where(np.array(LC[iPA][iEn]>=LbinsArr[iL]-0.125) & np.array(LC[iPA][iEn]<LbinsArr[iL+1]-0.125))[0]
                           MLTindex=np.where(np.array(MLTC[iPA][iEn]>=MLTbinsArr[iMLT]-0.25) & np.array(MLTC[iPA][iEn]<MLTbinsArr[iMLT+1]-0.25))[0]
                           # get the matches
                           matches=list(set(Lindex) & set(MLTindex))
                
                           PA_Flkp[iPA][iEn][iMLT][iL]+=list(np.array(FluxC[iPA][iEn])[matches])
                           
                           # this has not been screened for high or excessively low values
                           # so correction needs to applied when loading these values
                 # Check for month
  
        except:
                print 'error!'
        dt1=dt1+datetime.timedelta(days=1)
        yrCur=dt1.year
        print dt1
        if dt1.month != monthCur:
            # save this file
                os.chdir('/Users/loisks/Desktop/PA_check')
                if monthCur==12:
                    yrCur=yrCur-1
                for iEn in range(n1,n2):
                  for iPA in range(PAs):                      
                    pickling.hdf5_data_save(PA_Flkp[iPA][iEn],'Flux_sat=rbsp'+str(satellite[idir])+'_'+str(yrCur)+'_'+str(monthCur)+'_Kp='+str(kpBar)+'_En='+str(HopeEnergies[iEn])+'_PA='+str(PA_LABEL[iPA])+'_spe='+name_species[ispe] , 'CorrectedPAFluxes', Lbins,mltbins)
                monthCur=dt1.month
                print monthCur
            
                PA_Flkp=[ [ [  [ [] for x in xrange(Lbins)] for y in xrange(mltbins)] for z in  range(n1,n2)] for w in range(PAs)] # Energy x PA x MLT x L
                       


