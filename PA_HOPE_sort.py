# PA_HOPE_sort.py
#
# Sort by MLT, L-SHELL, and Pitch angle
# this will, unfortunately, take a while :( 
#
# LKS December 31st 2014. LAST CODE OF THE YEAR! LETS DO IT!
# SHOULD CHANGE THIS TO BE PA_HOPE_sort_MLAT.py if a re-run of all of this is necessary (which it probably will be)
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
os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/PA_Distributions')
#
# parameters
date1='20130201'
date2='20141001'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
Lbins=11
mltbins=48
LbinsArr=np.linspace(1.5, 4,nLbins)
MLTbinsArr=np.linspace(0.25, 23.75, nmlt_bins)
kpBar=3
# energy channels
n1=0
n2=16
PAs=11 # 11 pitch angles on HOPE
dir=['HOPE_A', 'HOPE_B']
name=['rbspa', 'rbspb']
satellite=['a','b']
species=[ 'H', 'He', 'O']
HopeEnergies=[.99,1.2,1.34,1.55,1.83,2.18,2.53,2.95,3.38,3.94,4.64,5.35,6.26,7.32,8.51,9.92]
name_species=['FPDU', 'FHEDU','FODU'] # gives the unidirectional flux with PA
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

                    
def weightedAvg(window,df,weight,DT):
      # window = time window
      # df = array to avg over
      # weight = array of weights
      # DT = starting date time
      cutoff=1440/window # number of windows in a day by minutes
      newDF=np.zeros(cutoff)
      weight=np.array(weight)
      for iTime in range( cutoff-1):
         sI=df.index.searchsorted(DT+iTime*datetime.timedelta(minutes=1))
         eI=df.index.searchsorted(DT+(iTime+1)*datetime.timedelta(minutes=1))
         # we have five minute window now get the weights
         totalC=1.0*np.nansum(weight[sI:eI])
         newDF[iTime]=np.nansum(np.array(df)[sI:eI]*weight[sI:eI])/totalC
      newDF[-1]=df[-1]
      timeArr=pd.date_range(DT,DT+datetime.timedelta(days=1), freq='1Min')[:-1]
      DF=pd.DataFrame(newDF, index=timeArr)
      # 1 minute weighted average
      return DF
#
#
# Need to set up that it saves every month
# Save by pitch angle
# save by s/c corrected energy
#    
# energy channels
for ispe in range(len(name_species)):
  for idir in range(len(dir)):
   # start overarching loop through days
   PA_Flkp=[ [ [  [ [] for x in xrange(Lbins)] for y in xrange(mltbins)] for z in range(PAs)] for w in range(n1,n2)] # Energy x PA x MLT x L
#   os.chdir('/Users/loisks/Desktop/PA_Check/CorrectedPAFluxes/')
#   check=glob.glob()
#   try:
#       check[0] # see if anything is there ,this will result in an index error if nothing is there
#       print species[ispe]+'_'+name[idir]+'_'+str(ien)
#   except(IndexError):
   dt1=dt0
   while dt1 < dt2:
        monthCur=dt1.month
        while True:
            try:
                date1=datetime.datetime.strftime(dt1, '%Y%m%d')
                os.chdir('/Users/loisks/Desktop/liemohn10/'+dir[idir])
                hope_file=name[idir]+'_rel02_ect-hope-PA-L3_'+date1+'*.cdf'
                ghope=glob.glob(hope_file)
                pyf=pycdf.CDF(ghope[0])
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+satellite[idir])
                f=glob.glob('*'+date1+'*')
                pyfe=pycdf.CDF(f[0])
                break
            except(IndexError):
                print dt1
                dt1=dt1+datetime.timedelta(days=1)
                pass
        # get the data quantities we need, separate by species 
                # now a 2 d array, chosen for correct PA
                
            # now resample and rebin           
            L=pyf['L_Ion'][...] # L Shell
            energy=pyf['HOPE_ENERGY_Ion'][...] # Hope Energy
            MLT=pyf['MLT_Ion'][...] # MLT
            epoch=pyf['Epoch_Ion'][...]    
            PA_LABEL=pyf['PITCH_ANGLE'][...][0] # right PA label
            
            #
            # resample by a minute
            rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
            rng=rt[::1]
            dfFlux={}
            dfCounts={}
            dfEnergy={}
            wMLT={}; wL={}; wEnergy={}; wFlux={}
            counter1=0
            kyotodf2=pd.DataFrame({'kp':nKP}, index=rng)
            kpdf=np.array(kyotodf2['kp'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
            highKp=np.where(kpdf>=kpBar)[0]
            # goes by pitch angle and initial energy channel
            for iPA in range(PAs):
                for iEn in range(n1,n2):
                    dfFlux[str(counter1)]=np.swapaxes(pyf[name_species[ispe]][...],2 ,0)[iEn][iPA]
                    dfCounts[str(counter1)]=np.swapaxes(pyf['Counts_'+spe[iSpe]][...],2,0)[iEn][iPA]
                    dfEnergy[str(counter1)]=np.swapaxes(energy, 1,0)[iEn]
                    #
                    # weighted averages
                    wMLT[str(counter1)]=weightedAvg(1, MLT, dfCounts[str(counter1)], dt1)
                    wL[str(counter1)]=weightedAvg(1,L,dfCounts[str(counter1)],dt1)
                    wEnergy[str(counter1)]=weightedAvg(1,dfEnergy[str(counter1)],dfCounts[str(counter1)],dt1)
                    wFlux[str(counter1)]=weightedAvg(1,dfFlux[str(counter1)],dfCounts[str(counter1)],dt1)
                    # nan the high Kps
                    wFlux[str(counter1)][highKp]=np.nan
                    #
                    # now nan high kp times
                    counter1+=1
            
            #
            # now get the EFW data and sort
            
            potential=-1*pyfe['Vavg'][...]
            pEpoch=pd.DatetimeIndex(pyfe['epoch'][...])
            df=pd.DataFrame(potential, index=pEpoch, columns=['potential'])
            Phi=np.array(df['potential'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
            Phi=np.absolute(Phi)
            # nan the times of significant s/c charging
            Phi[Phi>10]=np.nan
            counter2=0
            # now we need to rebin the fluxes based on s/c corrected energy
            FluxC=[ [ [] for i in range(n1,n2)] for j in range(PAs)]
            MLTC=[ [ [] for i in range(n1,n2)] for j in range(PAs)]
            LC=[ [ [] for i in range(n1,n2)] for j in range(PAs)]
            #
            # everything is on a 1440 scale now
            for iPA in range(PAs):
                for iEn in range(n1,n2):
                    wEnergy[str(counter2)]=np.array(wEnergy[str(counter2)])+np.array(Phi)
                    for iTime in len(wEnergy[str(counter2)]):
                        cEn=min(range(len(HopeEnergies)), key=lambda i: abs(HopeEnergies[i]-wEnergy[str(counter2)][iTime]))
                        FluxC[iPA][cEn].append(wFlux[str(counter2)][iTime])
                        MLTC[iPA][cEn].append(wMLT[str(counter2)][iTime])
                        LC[iPA][cEn].append(wL[str(counter2)][iTime])
                        # this should give an index error occasionally
 
                
             cflux[cflux > 1e12] = np.nan
             cflux[cflux <= 0] = np.nan
             for iPA in range(PAs):
              for iEn in range(n1,n2):
                 for iMLT in range(mltbins):
                   for iL in range(Lbins):
                       Lindex=np.where((np.array(LC[iPA][iEn]>=LbinsArr[iL]-0.125) & (LC[iPA][iEn]<LbinsArr[iL+1]-0.125))[0]
                       MLTindex=np.where((MLTC[iPA][iEn]>=MLTbinsArr[iMLT]-0.25) & (MLTC[iPA][iEn]<MLTbinsArr[iMLT+1]-0.25))[0]
                       # get the matches
                       matches=list(set(Lindex) & set(MLTindex))
    
                       PA_Flkp[iPA][iEn][imlt][iL]+=list(FluxC[iPA][cEn][matches])
                       # this has not been screened for high or excessively low values
                       # so correction needs to applied when loading these values
             # Check for month
             yrCur=dt1.year
             dt1=dt1+datetime.timedelta(days=1)
             if dt1.month != monthCur:
      # save this file
                os.chdir('/Users/loisks/Desktop/PA_check')
                for iEn in range(n1,n2):
                  for iPA in range(PAs):
                    pickling.hdf5_data_save(PA_Flkp[iPA][iEn],'Flux_sat=rbsp'+str(satellite[idir])+'_'+str(yrCur)+'_'+str(monthCur)+'_Kp='+str(kpBar)+'_En='+str(HopeEnergies[iEn])+'_PA='+str(PA_LABELS[iPA]) , 'CorrectedPAFluxes', Lbins, mltbins)
                    monthCur=dt1.month
                    print monthCur
            
             PA_Flkp=[ [ [  [ [] for x in xrange(Lbins)] for y in xrange(mltbins)] for z in range(PAs)] for w in range(n1,n2)] # Energy x PA x MLT x L
                       


