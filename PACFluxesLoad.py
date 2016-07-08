# PACFluxesLoad.py
#
# load in the PA corrected fluxes
# and plot them!
#
# LKS, November 2015, working towards that thesis!
#
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
import matplotlib.pyplot as plt
from matplotlib import cm
from dateutil.relativedelta import relativedelta
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
os.chdir('/Users/loisks/Desktop/Functions')
import pickling as pickling
import plots
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SummedPitchAngles')

#
date1='20130201'
date2='20150401'
#date2='20141001'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')

def diff_month(d1, d2):
    return (d1.year - d2.year)*12 + d1.month - d2.month
total_month=diff_month(dt2, dt0)
Lbins=11
mltbins=48
LbinsArr=np.linspace(1.5, 4.0,Lbins)
MLTbinsArr=np.linspace(0.25, 24.25,mltbins+1)
kpBarArr=[3]
# energy channels
n1=0
n2=16
PAs=11 # 11 pitch angles on HOPE
dir=['HOPE_A', 'HOPE_B']
name=['rbspa', 'rbspb']
satellite=['A','B']
species=[ 'P']#, 'He', 'O']
HopeEnergies=[.99,1.2,1.34,1.55,1.83,2.18,2.53,2.95,3.38,3.94,4.64,5.35,6.26,7.32,8.51,9.92]
PA_LABEL=[4.5,18.0,36.0,54.0,72.0,90.0,108.0,126.0,144.0,162.0,175.5]
PA_L=[0,9,27, 45, 63, 81, 99, 117, 135, 153, 171,180]
dt1=dt0

#
#
for iKP in range(len(kpBarArr)):
 os.chdir('CorrectedPAFluxes')
 kpBar=kpBarArr[iKP]
 PA_Flkp1=[ [ [  [ [] for x in range(Lbins)] for y in range(mltbins)] for z in range(n1,n2)] for w in range(PAs)]
 noPoints=[ [ [  [ 0 for x in range(Lbins)] for y in range(mltbins)] for z in range(n1,n2)] for w in range(PAs)]
 PA_Flkp3=[ [ [  [ [] for x in range(Lbins)] for y in range(mltbins)] for z in range(n1,n2)] for w in range(PAs)]
 mPA_Flkp1=[ [ [  [ [] for x in range(Lbins)] for y in range(mltbins)] for z in range(n1,n2)] for w in range(PAs)]
 mPA_Flkp3=[ [ [  [ [] for x in range(Lbins)] for y in range(mltbins)] for z in range(n1,n2)] for w in range(PAs)]
 meanPA_Flkp3=[ [ [  [ [] for x in range(Lbins)] for y in range(mltbins)] for z in range(n1,n2)] for w in range(PAs)]
 std_PA=[ [ [  [ [] for x in range(Lbins)] for y in range(mltbins)] for z in range(n1,n2)] for w in range(PAs)]

 
 dt1=dt0
 monthCount=1
 while dt1 < dt2:    
    monthCur=dt1.month
    yrCur=dt1.year
    print('Month Count: ' + str(monthCount) + ' / ' + str(total_month))
  #  if monthCur==12:
  #      yrCur+=1
    for iPA in range(len(PA_LABEL)):
        for iEn in range(len(HopeEnergies)):
         try:
            dataTemp1=pickling.hdf5_data_open('Flux_sat=rbspA_'+str(yrCur)+'_'+str(monthCur)+'_Kp='+str(kpBar)+'_En='+str(HopeEnergies[iEn])+'_PA='+str(PA_LABEL[iPA])+'_spe=FPDU.h5' , Lbins, mltbins)
            dataTemp2=pickling.hdf5_data_open('Flux_sat=rbspB_'+str(yrCur)+'_'+str(monthCur)+'_Kp='+str(kpBar)+'_En='+str(HopeEnergies[iEn])+'_PA='+str(PA_LABEL[iPA])+'_spe=FPDU.h5' , Lbins, mltbins)
            for iMLT in range(mltbins):
                for iL in range(Lbins):
                    T1a=np.array(dataTemp1[iMLT][iL])
                    PA_Flkp3[iPA][iEn][iMLT][iL]+=list(T1a[~np.isnan(T1a)])
                    T2=np.array(dataTemp2[iMLT][iL])
                    PA_Flkp3[iPA][iEn][iMLT][iL]+=list(T2[~np.isnan(T2)])
         except:
          print ('No Data!')
                    #if len(list(T2[~np.isnan(T2)])) < len(list(T2[~np.isnan(T2)])):
                    #        print "BAD"
                    #PA_Flkp[iPA][iEn][iMLT][iL]+=list(dataTemp[iMLT][iL])#+list(dataTemp2[iMLT][iL])
    dt1=dt1+relativedelta(months=+1)
    monthCount+=1
#
# now get median
 print('getting medians and standard deviation') 
 for iPA in range(len(PA_LABEL)):
    for iEn in range(len(HopeEnergies)):
        for iMLT in range(mltbins):
            for iL in range(Lbins):
                # remove the zeros                
                PA_Flkp3[iPA][iEn][iMLT][iL]=np.array(PA_Flkp3[iPA][iEn][iMLT][iL])              
                PA_Flkp3[iPA][iEn][iMLT][iL][np.where(PA_Flkp3[iPA][iEn][iMLT][iL]==0)[0]]=np.nan
                mPA_Flkp3[iPA][iEn][iMLT][iL]=np.nanmedian(PA_Flkp3[iPA][iEn][iMLT][iL])
                meanPA_Flkp3[iPA][iEn][iMLT][iL]=np.nanmean(PA_Flkp3[iPA][iEn][iMLT][iL])
                std_PA[iPA][iEn][iMLT][iL]=np.nanstd(PA_Flkp3[iPA][iEn][iMLT][iL])
                noPoints[iPA][iEn][iMLT][iL]=len(PA_Flkp3[iPA][iEn][iMLT][iL])
                #if (np.isnan(mPA_Flkp3[iPA][iEn][iMLT][iL]) == True) or (mPA_Flkp3[iPA][iEn][iMLT][iL]==0):
                     #mPA_Flkp3[iPA][iEn][iMLT][iL]=mPA_Flkp1[iPA][iEn][iMLT][iL]
                    

#then plot
 lScales=[8,8,8,8,7,7,7,7,6,6,5,5,5,5,5,5]
 hScales=[11,11,11,11,10,10,10,10,9,9,8,8,8,8,8,8]

 os.chdir('/Users/loisks/Documents/ResearchProjects/SummedPitchAngles')
 for iPA in range(len(PA_LABEL)):
    for iEn in range(len(HopeEnergies)):
        plots.fancy_plot(mPA_Flkp3[iPA][iEn], mltbins, Lbins, 6, 10, 'Differential Number Flux', 'H$^{+}$', str(kpBar), '_noZeros_PA='+str(PA_LABEL[iPA])+'_En='+str(HopeEnergies[iEn]), 'PACorrectedPlots', 'log', 'jet')
        plots.fancy_plot(noPoints[iPA][iEn], mltbins, Lbins, 0, 4, 'Number of Points', 'H$^{+}$', str(kpBar), '_noZeros_PA='+str(PA_LABEL[iPA])+'_En='+str(HopeEnergies[iEn]), 'NumPoints', 'log', 'viridis')
        plots.fancy_plot(np.array(std_PA[iPA][iEn])/np.array(meanPA_Flkp3[iPA][iEn]), mltbins, Lbins, 0, 4, '$\sigma_{D}$ / $\mu$', 'H$^{+}$', str(kpBar), 'STD_noZeros_PA='+str(PA_LABEL[iPA])+'_En='+str(HopeEnergies[iEn]), 'STD_PACorrectedPlots', '', 'Purples')
        #plots.fancy_plot(mPA_Flkp1[iPA][iEn], mltbins, Lbins, 6, 10, 'Differential Number Flux', 'H$^{+}$', str(kpBar), '_noZeros_PA='+str(PA_LABEL[iPA])+'_En='+str(HopeEnergies[iEn]), 'PACorrectedPlots', 'log', 'jet')
 stop

 for iEn in range(len(HopeEnergies)):
    lowScale=lScales[iEn]
    highScale=hScales[iEn]
    for iL in range(Lbins):
                L_label=str(0.25*iL+1.5)
                #mm1=np.swapaxes(mPA_Flkp1, 0, 1,)[iEn] # gives PA x MLT x L
                mm3=np.swapaxes(mPA_Flkp3, 0, 1)[iEn]
                mmstd=np.swapaxes(std_PA,0,1)[iEn]
                mean3=np.swapaxes(meanPA_Flkp3, 0, 1)[iEn]
                #data_L1 = np.swapaxes(mm1, 2, 0)[iL] # gives MLT x PA
                data_L3 = np.swapaxes(mm3, 2, 0)[iL]
                dd=np.swapaxes(mean3, 2, 0)[iL]
                #dataSTD= 10**(np.array(np.swapaxes(mmstd,2,0)[iL])+np.log10(data_L3))/data_L3
                dataSTD=np.array(np.swapaxes(mmstd,2,0)[iL])/dd
                xdata_size = np.linspace(0, 24, mltbins)
                ydata_size = PA_LABEL
                EN_LABEL=HopeEnergies[iEn]
                #plots.box_plot(data_L1, xdata_size, PA_L, 0, 24,0, 180, 'MLT (hours)', 'PA',lowScale, highScale, 'Differential # Flux', species[0], 'PA_vs_MLT' ,'corrected_'+ species[0]+'noZero_kp='+str(1)+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA')
                #plots.box_plot(data_L1, xdata_size, PA_L, 0, 6,0, 180, 'MLT (hours)', 'PA', lowScale, highScale, 'Differential # Flux', species[0], 'PA_vs_MLT_0_6' , 'corrected_'+species[0]+'noZero_kp='+str(1)+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA')
                plots.box_plot(data_L3, xdata_size, PA_L, 0, 24,0, 180, 'MLT (hours)', 'PA',lowScale, highScale, 'Differential # Flux', species[0], 'SC_PA_vs_MLT' ,'corrected_'+ species[0]+'noZero_kp='+str(3)+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA')
                plots.box_plot(data_L3, xdata_size, PA_L, 0, 6,0, 180, 'MLT (hours)', 'PA', lowScale, highScale, 'Differential # Flux', species[0], 'SC_PA_vs_MLT_0_6' , 'corrected_'
+species[0]+'noZero_kp='+str(3)+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA')
                plots.linear_box_plot(dataSTD, xdata_size, PA_L, 0, 24,0, 180, 'MLT (hours)', 'PA',0, 3, '$\sigma_{D}$ / $\mu$', species[0], 'STD_SC_PA_vs_MLT' ,'STD_corrected_'+ species[0]+'noZero_kp='+str(3)+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA')#, colors='Purples')

for iL in range(Lbins):
# for imlt in range(mltbins):
   mm3=np.swapaxes(mPA_Flkp3, 0, 3,)[iL] # PA x MLT x EN
   #mm4=np.swapaxes(mm3, 0,1)[imlt] # this leaves PA x EN
   xdata_size = PA_LABEL
   ydata_size = HopeEnergies
   colorbar_max=10
   colorbar_min=6
   Llabel=str(LbinsArr[iL])
   #MLTlabel=str(MLTbinsArr[imlt])
   plt.register_cmap(name='viridis', cmap=cmaps.viridis)
   plt.set_cmap(cmaps.viridis)
   #
   # code for special log box plot
   from numpy import ma
   #cmap=plt.cm.viridis
   # intitialize the figure
   fig=plt.figure()
   MLTs=['2', '4', '6', '8', '10', '12', '14', '16', '18','20','22','0']
   nMLT=[4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 0]
   
   fig, ((ax1, ax2,ax3, ax4), (ax5, ax6,ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
   fig.subplots_adjust(hspace=.25)
   combo=[ax1, ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]
   for imlt in range(len(MLTs)):
     ax=combo[imlt]
     mm4=np.swapaxes(mm3, 0,1)[nMLT[imlt]]
     #ax.set_xlabel('PA [degrees]', fontsize=30, fontweight='bold')
     #ax.set_ylabel('Energy [eV]', fontsize=30, fontweight='bold')
     plt.subplots_adjust(right=0.8, top=0.9, bottom=0.15)
     # mask the nans
     # log mods here 
     data=np.log10(np.array(mm4))
     datah_m=ma.masked_invalid(data)
     X,Y=np.meshgrid(xdata_size, ydata_size)
     col=ax.pcolormesh(X,Y,datah_m,cmap='viridis', vmax=colorbar_max, vmin=colorbar_min)
     ax.set_xlim(5,175)
     ax.set_ylim(1.5,10)
     ax.set_yscale('log')
     ax.set_title('MLT = ' + str(MLTs[imlt]), fontsize=22, fontweight='bold')
     ax.tick_params(axis='both', which='major', labelsize=25)     
     ax.set_xticks([5, 90, 175])
     ax.set_yticks([1.5, 10])
     ax.set_yticklabels(['1.5', '10'])
     ax.set_xticklabels(['5', '90', '175'], fontsize=22, fontweight='bold')
   cbaxes = fig.add_axes([0.81, 0.15, 0.03, 0.75]) 
   cb =fig.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(colorbar_min, colorbar_max+1) ) 
   cb.set_label('cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$', fontsize=30, fontweight='bold')
   cb.ax.tick_params(labelsize=30) 
   #  ax.set_title('MLT = ' + str(MLTlabel), fontsize=25)
   subdir_name='Energy_PA_Box'
   plt.draw()
   if not os.path.exists(subdir_name):
       os.umask(0) # unmask if necessary
       os.makedirs(subdir_name) 
   os.chdir(subdir_name)#
   fig.set_size_inches(18,9)
   plt.savefig('EnergyPAGrid_L='+Llabel+'.png')
   plt.close(fig)
   os.chdir('..')
#
n1a=3
n2a=12
cmap=plt.cm.jet
# create line plots at one L Shell for all energies and MLTS for PA = 90
for iL in range(Lbins):
    dataL=np.swapaxes(np.swapaxes(mPA_Flkp3, 1, 3)[5][iL],1,0) # L = 2 , PA = 90 l eaves  MLT x ien
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.subplots_adjust(right=0.75, top=0.85, bottom=0.15)
    ax.set_xlabel('MLT', fontsize=25, fontweight='bold')
    ax.set_ylabel(' Normalized # Flux', fontsize=25, fontweight='bold')
    ax.set_yscale('log')
    ax.set_xlim(0, 24)
    ax.set_ylim(1e-3,1e0)
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
    plt.rc('font', **font)
    # data and plotting loop
    plt.gca().set_color_cycle([cmap(i) for i in np.linspace(0, 0.99, n2a-n1a)])
    labels=[]
    for iy in range(n1a,n2a):
        ax.plot(range(mltbins), dataL[iy]/np.nanmax(dataL[iy]), lw=3, marker='x')
        labels.append( str(HopeEnergies[iy]) + ' eV')
    plt.legend(labels, bbox_to_anchor=[1.37, 1], fontsize=25)
    subdir_name='SC_lineplots'
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name) 
    os.chdir(subdir_name)#
    fig.set_size_inches(13,9)
    plt.savefig('lineplot_L='+str(iL)+'.png')
    plt.close(fig)
    os.chdir('..')

for iPA in range(len(PA_LABEL)):
    dataL=np.swapaxes(np.swapaxes(mPA_Flkp3, 1, 3)[iPA][4],1,0) # L = 2 , PA = 90 l eaves  MLT x ien
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.subplots_adjust(right=0.75, top=0.85, bottom=0.15)
    ax.set_xlabel('MLT', fontsize=25, fontweight='bold')
    ax.set_ylabel(' Normalized # Flux', fontsize=25, fontweight='bold')
    ax.set_yscale('log')
    ax.set_xlim(0, 24)
    ax.set_ylim(1e-3,1)
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
    plt.rc('font', **font)
    # data and plotting loop
    plt.gca().set_color_cycle([cmap(i) for i in np.linspace(0, 0.99, n2a-n1a)])
    labels=[]
    for iy in range(n1a,n2a):
        # figure out why this is not showing 1
        alpha=np.array( dataL[iy])/np.nanmax(np.array(dataL[iy]))
        ax.plot(range(mltbins),alpha, lw=3, marker='x')
        labels.append( str(HopeEnergies[iy]) + ' eV')
    plt.legend(labels, bbox_to_anchor=[1.37, 1], fontsize=25)
    subdir_name='SC_lineplots'
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name) 
    os.chdir(subdir_name)#
    fig.set_size_inches(13,9)
    plt.savefig('lineplot_PA='+str(iPA)+'.png')
    plt.close(fig)
    os.chdir('..')
