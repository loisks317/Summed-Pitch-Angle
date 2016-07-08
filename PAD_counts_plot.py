# PAD_counts_plot.py
#
# takes the total counts in each energ bin and plots the number
# 
# LKS November 2015. Day before Thanksgiving. At least it's not Thanksgiving.
# Did end up finishing this on Thanksgiving
#
import numpy as np
from spacepy import pycdf
import glob
import os
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
import pickle 
import matplotlib.pyplot as plt
import h5py
from scipy import interp
from dateutil.relativedelta import relativedelta
import os
import numpy as np
from matplotlib import pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
from numpy import ma
# INPUT PARAMETERS
#
# start date
date1='20130201'
#
# end date
date2='20150401'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
    # gives the number of months in time range for easy file output
    return (d1.year - d2.year)*12 + d1.month - d2.month
total_month=diff_month(dt2, dt0)
#
# parameters 
Lbins=11
mltbins=48
Lmax=4
Lmin=1.5
# energy channels
n1=0 # 0.99 eV
n2=16 # 9.9 eV
PAs=11 # 11 pitch angles on HOPE
Energies=[.99,1.2,1.34,1.55,1.83,2.18,2.53,2.95,3.38,3.94,4.64,5.35,6.26,7.32,8.51,9.92]
name=['rbspa', 'rbspb']
#
# bins to sum
window=10
dataArr={}
for ien in range(len(Energies)):
  dt1=dt0
  data=[ [[ 0 for i in range(Lbins)] for j in range(mltbins)] for k in range(12)]
  os.chdir('Window'+str(window))    
  for imonth in range(total_month):
    cur_date=str(dt1.month)+'_'+str(dt1.year)

    
    for iCount in range(12):
        combo=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
        temp=pickle.load(open(cur_date+'_sc_counts='+str(iCount)+'_H_'+str(ien)+'_rbspa_'+str(Energies[ien])+'_countsHisto.p', 'rb'))
        tempb=pickle.load(open(cur_date+'_sc_counts='+str(iCount)+'_H_'+str(ien)+'_rbspb_'+str(Energies[ien])+'_countsHisto.p', 'rb'))
        for imlt in range(mltbins):
            for iL in range(Lbins):
                combo[imlt][iL]=temp[imlt][iL]+tempb[imlt][iL]
                data[iCount][imlt][iL]+=combo[imlt][iL]
      
        
    dt1=dt1+relativedelta(months=1)
  os.chdir('..')  
  dataPlot=np.zeros(12)
  # now sum all MLTs and take the L = 2.5 one
  for iC in range(12):
    dataPlot[iC]=np.sum(data[iC], axis=0)[2] # L = 2.5
    #
    # plot this for each energy channel
  fig=plt.figure()
  ax=fig.add_subplot(111)
  fig.subplots_adjust(bottom=0.15, left=0.16, right=0.79)
  font = {'family' : 'normal',
              'weight' : 'bold',
              'size'   : 22}
  plt.rc('font', **font)
      # combine the bad counts and uncategorized
  
      #ax.set_title(plot_title + ' at L=' + str(L_cur))
  ax.set_ylabel('Counts Across all MLTs', fontweight='bold')
  ax.set_xlabel('Number of Invalid Indices', fontweight='bold')
  ax.set_xlim(-0.6, 11.6)
  ax.set_ylim(10,10000)
  ax.set_yscale('log')
  r1=ax.bar(range(12), dataPlot, color='DarkViolet',align='center')
  params = {'legend.fontsize': 15}
  plt.rcParams.update(params)
  plt.grid()
  dataArr[ien]=dataPlot


  
  subdir_name='PACountsBarPlots'
  if not os.path.exists(subdir_name):
   os.umask(0) # unmask if necessary
   os.makedirs(subdir_name, 0777) 
  os.chdir('PACountsBarPlots')#
  plt.savefig('PACounts_L=2.5_MLTsummed_n='+str(ien)+'.png')
  plt.close()
  os.chdir('..')

# now triple plot
fig=plt.figure()
fig.subplots_adjust(hspace=.4)
#plt.figure(figsize=(20,10))
ax1=fig.add_subplot(311)
fig.subplots_adjust(bottom=0.15, left=0.16, right=0.9)
font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
plt.rc('font', **font)

#ax.set_ylabel('Counts Across all MLTs', fontweight='bold')
#ax.set_xlabel('Number of Bad Indices', fontweight='bold')
ax1.set_xlim(-0.6, 11.6)
ax1.set_yscale('log')
ax1.set_ylim(1e-3,1e0)
ax1.set_title('1.5 eV', fontweight='bold')
d3=np.array(dataArr[3])
r1=ax1.bar(range(12), d3/np.nanmax(d3), color='DarkViolet',align='center')
ax1.set_xticklabels([])
plt.rcParams.update(params)
plt.grid()
dataArr[ien]=dataPlot
ax2=fig.add_subplot(312)
ax2.set_ylabel('Normalized Counts \n Summed Across all MLTs', fontweight='bold', fontsize=25)
#ax.set_xlabel('Number of Bad Indices', fontweight='bold')
ax2.set_xlim(-0.6, 11.6)
ax2.set_yscale('log')
ax2.set_title('3.0 eV', fontweight='bold')
dd7=np.array(dataArr[7])
ax2.set_ylim(1e-3,1e0)
r1=ax2.bar(range(12), dd7/np.nanmax(dd7), color='DarkViolet',align='center')
#ax2.set_ylim(-1e3,1e1)
ax2.set_xticklabels([])
plt.rcParams.update(params)
plt.grid()
dataArr[ien]=dataPlot
ax3=fig.add_subplot(313)
#ax3.set_ylabel('Counts Across all MLTs', fontweight='bold')
ax3.set_xlabel('Number of Invalid Indices', fontweight='bold', fontsize=25)
ax3.set_xlim(-0.6, 11.6)
ax3.set_ylim(1e-3,1e0)
ax3.set_yscale('log')
ax3.set_title('5.3 eV', fontweight='bold')
d11=np.array(dataArr[11])
r1=ax3.bar(range(12), d11/np.nanmax(d11), color='DarkViolet',align='center')
plt.xticks(range(12))
plt.rcParams.update(params)
plt.grid()
dataArr[ien]=dataPlot
fig.set_size_inches(12, 9)
subdir_name='PACountsBarPlots'
if not os.path.exists(subdir_name):
 os.umask(0) # unmask if necessary
 os.makedirs(subdir_name, 0777) 
os.chdir('PACountsBarPlots')#
plt.savefig('triple_PACounts_L=2.5_MLTsummed_n='+str(ien)+'.png')
plt.close()
os.chdir('..')
