# PA_sorted_plot.py
#
# take the output of PA_HOPE_sort.py (should be in this folder) and assimilate
# and plot the results for large statistical study
#
# LKS January 2015 [5 MONTHS TILL QUAL]
#
# imports
import numpy as np
from spacepy import pycdf
import glob
import os
import pickle
import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
import matplotlib.collections as collections
import matplotlib.dates as dt
from matplotlib.colors import LogNorm
os.chdir('/Users/loisks/Documents/Functions')
import plots
import pickling
os.chdir('/Users/loisks/Desktop/PA_Distributions')
#
# globals
ntotal=16
PAs=11
Lbins=26
mltbins=49
Energy_labels=[0.99, 1.2,1.33, 1.5, 1.8, 2.2, 2.53, 2.95, 3.38, 3.9, 4.6, 5.34, 6.26, 7.31, 8.51, 9.9]
with open('PA_LABELS', 'rb') as f0:
    PA_labels=pickle.load(f0)
spe=['H', 'He', 'O']
for ispe in range(len(spe)):
    for n1 in range(ntotal):
        all_PA={}
# d = hdf5_data_open(filename)
# some_arr = [ [] of stuff of stuff ] 
# count=count+1
# for imlt in range(mltbins):
#  for iL in range(Lbins):
#     some_arr[imlt][iL]=d[str(count)]
#     count=count+1
     #   for iPA in range(len(PA_labels)):
        meanf=[[ [x for x in xrange(len(PA_labels))] for x in xrange(Lbins)] for x in xrange(mltbins)]
        no_points=[[ [x for x in xrange(len(PA_labels))] for x in xrange(Lbins)] for x in xrange(mltbins)]
        for iPA in range(len(PA_labels)):
#
# read in the data files
            print PA_labels[iPA]
            os.chdir('PA_sorted_flux')
            filename1=glob.glob(spe[ispe]+'_F_lkp_out_rbspa_'+str(n1)+'*'+str(PA_labels[iPA])+'.h5')
            filename2=glob.glob(spe[ispe]+'_F_lkp_out_rbspb_'+str(n1)+'*'+str(PA_labels[iPA])+'.h5')
            d=pickling.hdf5_data_open(filename1[0]) 
            d2=pickling.hdf5_data_open(filename2[0]) 
            file=[[ [] for x in xrange(Lbins)] for x in xrange(mltbins)]
            file2=[[ [] for x in xrange(Lbins)] for x in xrange(mltbins)]
            cfile=[[ [] for x in xrange(Lbins)] for x in xrange(mltbins)]
            count=0
            for imlt in range(mltbins):
                for iL in range(Lbins):
                    file[imlt][iL]=d[str(count)]    
                    file2[imlt][iL]=d2[str(count)]
# combine A and B
                    cfile[imlt][iL]=np.array(list(file[imlt][iL])+list(file2[imlt][iL]))
                    count=count+1
#
# okay now exclude anamolies 
            for imlt in range(mltbins):
                for iL in range(Lbins):
                    try:
                        cfile[imlt][iL]=cfile[imlt][iL][~np.isnan(cfile[imlt][iL])]
                        cfile[imlt][iL]=cfile[imlt][iL][cfile[imlt][iL] > 0]
                        cfile[imlt][iL]=cfile[imlt][iL][cfile[imlt][iL] < 1e14]
                        no_points[imlt][iL][iPA]=len(cfile[imlt][iL])
                        meanf[imlt][iL][iPA]=cfile[imlt][iL].mean() # screened for nans, zeroes, and high fluxes now I think
                    except(TypeError, AttributeError):
                        meanf[imlt][iL][iPA]=np.nan
                        no_points[imlt][iL][iPA]=0
            os.chdir('..')
#
# Now plot
 #           PA_LABEL=PA_labels[iPA]
        EN_LABEL=Energy_labels[n1]
        # now do an L-Shell loop
        for iL in range(5, 17):
            L_label=str(0.25*iL)
            data_L = np.swapaxes(meanf, 1, 0)[iL] # gives MLT x PA 
            nop=np.swapaxes(no_points, 1, 0)[iL]
            xdata_size = np.linspace(0, 24, mltbins)
            ydata_size = PA_labels
            plt_title = 'Energy = ' + str(EN_LABEL) + ' eV, L-Shell = ' + str(L_label)
            plots.box_plot(data_L, xdata_size, ydata_size, 0, 24,0, 180, 'MLT (hours)', 'PA', plt_title, 'Differential # Flux', spe[ispe], 'PA_vs_MLT' , spe[ispe]+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA')
            plots.box_plot(nop, xdata_size, ydata_size, 0, 24,0, 180, 'MLT (hours)', 'PA', plt_title, 'Number of Points', spe[ispe], 'PA_vs_MLT' , spe[ispe]+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA_npoints')
# zoom in on post midnight region
            plots.box_plot(data_L, xdata_size, ydata_size, 0, 6,0, 180, 'MLT (hours)', 'PA', plt_title, 'Differential # Flux', spe[ispe], 'PA_vs_MLT_0_6' , spe[ispe]+'_energy='+str(EN_LABEL)+'_Lshell='+L_label+'_MLT_vs_PA')


        for imlt in range(1,mltbins):
            mlt_label=str(0.5*imlt)
            data_mlt = meanf[imlt] # gives MLT x PA 
            nop=no_points[imlt]
            xdata_size = np.linspace(0,6.25, Lbins)
            ydata_size = PA_labels
            plt_title = 'Energy = ' + str(EN_LABEL) + ' eV, MLT = ' + str(mlt_label)
            plots.box_plot(data_mlt, xdata_size, ydata_size, 1.5, 6.25,0, 180, 'L-Shell', 'PA', plt_title, 'Differential # Flux', spe[ispe], 'PA_vs_LShell' , spe[ispe]+'_energy='+str(EN_LABEL)+'_MLT='+mlt_label+'_LShell_vs_PA')
            plots.box_plot(nop, xdata_size, ydata_size, 1.5, 6.25,0, 180, 'Lshell ', 'PA', plt_title, 'Number of Points', spe[ispe], 'PA_vs_Lshell' , spe[ispe]+'_energy='+str(EN_LABEL)+'_MLT='+L_label+'_Lshell_vs_PA_npoints')
     #   os.chdir('..')



           # plots.fancy_plot(meanf, mltbins, Lbins, 5, 9, 'Differential # Flux', spe[ispe], 'low', 'PA_stats_'+str(n1)+'_'+str(PA_labels[iPA]), 'PA_stats', 'log','PA = ' + str(PA_LABEL) + ', Energy = '+ str(EN_LABEL))
           # all_PA[iPA]=meanf




