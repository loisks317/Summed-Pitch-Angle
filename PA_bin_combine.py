# PA_Bin_combine.py
#
# make sure things work out count wise so we have good 
# statistics for our binning
#
# LKS March 2015, edited June 2015
# modification of PA_Bin_sort.py, making sure everything is correct
# and taking out the mu space correction
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
import matplotlib.pyplot as plt
import h5py
os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/PA_check')
from scipy import interp
from dateutil.relativedelta import relativedelta
#
# parameters
#date1='20130201'
#date2='20141002'
#
# just one month for now
date1='20130201'
date2='20140801'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
    return (d1.year - d2.year)*12 + d1.month - d2.month
total_month=diff_month(dt2, dt0)
Lbins=26
mltbins=49
# energy channels
n1=10
n2=13
PAs=11 # 11 pitch angles on HOPE
#
# bins to sum
window=10
PA_labels=[4.5, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0, 175.5]
mid_points=np.array(  [13.5, 45, 90, 135, 166.5])
dir=['HOPE_A', 'HOPE_B']
name=['rbspa', 'rbspb']
species=['H']#, 'He', 'O']
name_species=[ 'FPDU']#, 'FHEDU','FODU'] # gives the unidirectional flux with PA

os.chdir('/Users/loisks/Desktop/PA_check/PA_Combine_Bins/')
if not os.path.exists('Window='+str(window)):
        os.umask(0) # unmask if necessary
        os.makedirs('Window='+str(window), 0777) 
os.chdir('..')
#
# energy channels
for ispe in range(len(name_species)):
 for ien in range(n1,n2):
  for idir in range(len(dir)):
   dt1=dt0
   for imonths in range(total_month):
    # get the times right
    cur_date=str(dt1.month)+'_'+str(dt1.year)
    dt2=dt1+relativedelta(months=1)
    # start overarching loop through days
    cEn=[] # for the total energy average at the end
    delta=[]
# 
# DEFINE ARRAYS
    inverse_butterfly=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    source_cone=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    one_cone=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    butterfly=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    isotropic=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    loss_cone=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    bad_counts=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    uncat=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    spacecraft_charging=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    missing_fluxes=[[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)] for x in xrange(len(PA_labels))]
    # five is for the number of categories we sort the HOPE PAs by
    cats=[inverse_butterfly, source_cone, one_cone, butterfly, isotropic, loss_cone, uncat, bad_counts, spacecraft_charging]
    cat_labels=[ 'inverse_butterfly', 'source_cone', 'one_cone', 'butterfly', 'isotropic', 'loss_cone', 'uncategorized','bad_counts', 'spacecraft_charging']
    os.chdir('/Users/loisks/Desktop/PA_check/PA_Combine_Bins/')
    check=glob.glob(cur_date+'_'+species[ispe]+'_binned_fluxes_combine_'+ name[idir]+'_'+str(ien)+'*')
    try:
        check[0] # see if anything is there ,this will result in an index error if nothing is there
        print cur_date+'_'+species[ispe]+'_'+name[idir]+'_'+str(ien)
        dt1=dt2
    except(IndexError):
     while dt1 < dt2:
        while True:
            try:
                date1=datetime.datetime.strftime(dt1, '%Y%m%d')
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/'+dir[idir])
                hope_file=name[idir]+'_rel02_ect-hope-PA-L3_'+date1+'*.cdf'
                ghope=glob.glob(hope_file)
                pyf=pycdf.CDF(ghope[0])
                break
            except(IndexError):
                print dt1
                dt1=dt1+datetime.timedelta(days=1)
                pass
        # get the data quantities we need, separate by species 
                # now a 2 d array, chosen for correct PA
        L=pyf['L_Ion'][...] # L Shell
        energy=pyf['HOPE_ENERGY_Ion'][...] # Hope Energy
        MLT=pyf['MLT_Ion'][...] # MLT
        FLUX=pyf[name_species[ispe]][...]
        COUNTS=pyf['Counts_P'][...] # time x 11 x 72
        epoch=pyf['Epoch_Ion'][...]    
        PA_LABEL=pyf['PITCH_ANGLE'][...] # right PA label 
        os.chdir('/Users/loisks/Desktop/PA_check')
        kyoto=spk.fetch('kp',dt1, dt1)
        day=datetime.datetime.strftime(dt1,'%d') # need to get day of month
        kp_arr=kyoto['kp'][(int(day)-1)*8: int(day)*8]
        hours=np.array([0, 3, 6, 9, 12, 15, 18, 21])
        #
        temp_en=[]
        count=0
        for itime in range(len(epoch)/window): 
            while True: # everything ends in a break
               iLArr=np.array(L[itime*window:(window*itime+window)])/0.25
               iL=int(np.round(np.nanmean(iLArr)))
               if (iL >= 26): # if this is greater than zero
                   break # this means that we're at apogee which is beyond our area 
                   # of interest
               iMltArr=np.array(MLT[itime*window:(window*itime+window)])/0.5
               imlt=int(np.round(np.nanmean(iMltArr)))
               cengArr=np.swapaxes(np.array(energy[itime*window:(window*itime+window)]), 1, 0)[ien]# 10 array
               if len(np.where(cengArr > 20)[0]) > 0:
                   # perigee mode 
                   break
               cengArrMean=np.nanmean(cengArr)
               # Kp!
               middleTime=epoch[window*itime+(window/2)]
               hourIndex=np.where(hours <= middleTime.hour)[0][-1]
               if kp_arr[hourIndex] > 3:
                   break # high activity
               temp_en.append(cengArrMean)
               cfluxArr=np.array(np.swapaxes(FLUX[(window*itime):(itime*window+window)], 2, 0)[ien])
               cfluxArr[cfluxArr > 1e12]=np.nan # invalid stuff
               cfluxArr[cfluxArr < 0]=np.nan
               # PA x Time for both of these... I think
               cCountsArr=np.array(np.swapaxes(COUNTS[(window*itime):(itime*window+window)], 2, 0)[ien])
               # 
               # now sum the counts for each PA and give summary statistics and do weighted averages
               # of the flux based on the counts
               summedCounts=np.sum(cCountsArr, axis=1)
               notValidCounts=np.where(summedCounts < 10)[0] # define 10 as the cut off point
               # 
               # do a weighted average of the fluxex
               weightedMultiplication = cfluxArr*cCountsArr
               # now sum along the time axis
               summedWeights=np.nansum(weightedMultiplication, axis=1)
               weightedFlux=summedWeights/summedCounts
               weightedFlux[notValidCounts]=np.nan
               # probably rarely will have a full distribution now
               # due to counting issues
               Low_indices=np.where(np.isnan(weightedFlux)==True)[0]
               # change this to must have 5 good counts
               if len(Low_indices) >= 6:
                  bad_counts[imlt][iL]+=1
                  if (len(Low_indices) == 11):
                      spacecraft_charging[imlt][iL]+=1
                    #  print "sc charging all zeros at " + str(epoch[itime+window/2])
                  break # PAD is mostly zeros, just get rid of it
               # do not replace here 
               #try:
               #  weightedFlux[Low_indices]=np.min(weightedFlux[weightedFlux > 0]) ### NEW LOW VALUES  
               #except(ValueError):
               #  break
               #
               # now categorize where there are missing counts
               missing_indexes=np.where(np.isnan(weightedFlux)==True)[0]
               for im in (missing_indexes):
                  missing_fluxes[im][imlt][iL]+=1
                  # MISSING END1
                  if (weightedFlux[0]==np.nan) & (weightedFlux[1]==np.nan):
                      missingCats[0][imlt][iL]=missingCats[0][imlt][iL]+1
                  # MISSING INTR1
                  if (weightedFlux[2]==np.nan) & (weightedFlux[3]==np.nan):
                      missingCats[1][imlt][iL]=missingCats[1][imlt][iL]+1
                  # MISSING MIDDLE
                  if (weightedFlux[4]==np.nan) & (weightedFlux[5]==np.nan) & (weightedFlux[6]==np.nan):
                      missingCats[2][imlt][iL]=missingCats[2][imlt][iL]+1
                  # MISSING INTR2
                  if (weightedFlux[7]==np.nan) & (weightedFlux[8]==np.nan):
                      missingCats[3][imlt][iL]=missingCats[3][imlt][iL]+1
                  # MISSING END1
                  if (weightedFlux[9]==np.nan) & (weightedFlux[10]==np.nan):
                      missingCats[4][imlt][iL]=missingCats[4][imlt][iL]+1
                  break
               # normalize the flux by mean
               nfactor=np.nanmean(weightedFlux)
               nflux=weightedFlux/nfactor # now it's normalized
               zeros=len(np.where(nflux==0)[0])
               nans=len(np.where(np.isnan(nflux)==True)[0])
               if (zeros + nans) > 5:
                  # print "spacecrafting charging event"
                   spacecraft_charging[imlt][iL]+=1
                   print "sc charging nans and zeros at " + str(epoch[itime+window/2])
                   break

               fit=interp(np.linspace(0, 180, 181), np.array(PA_labels)[~np.isnan(nflux)], nflux[~np.isnan(nflux)])
               # okay try to define quartiles here
               # sort the array to exlude the highest and lowest points
               oflux = np.sort(nflux[~np.isnan(nflux)])
               min_value=oflux[1] # exclude outlier 
               peak_value=oflux[-2] # exclude outlier
               isoCheck=peak_value/min_value
            #   log_range = np.logspace(-2, 2, 1000)
            #   peak_a=np.where( log_range > peak_value)[0][0]
            #   min_a=np.where( log_range> min_value)[0][0]
               end1=np.nanmean(nflux[0:2])
               sec1=np.nanmean(nflux[2:4])
               middle=np.nanmean(nflux[4:7])
               sec2=np.nanmean(nflux[7:9])
               end2=np.nanmean(nflux[9:11])
               # 
               # sort thigns
               if isoCheck < 1.2: # if the maximum is less than 1.2x the minimum
                    isotropic[imlt][iL]=isotropic[imlt][iL]+1
                    plot_title='Isotropic'
                    count=count+1 
               elif (end1 < sec1) and (end2 < sec2) and (sec1 > middle) and (sec2 > middle):
                butterfly[imlt][iL]=butterfly[imlt][iL]+1
                plot_title='Butterfly'
                count= count+1
              # Inverse Butterfly
               elif (end1 > sec1) and (end2 > sec2) and (sec1 < middle) and (sec2 < middle):
                  inverse_butterfly[imlt][iL]=inverse_butterfly[imlt][iL]+1
                  plot_title='Inverse Butterfly'
                  count= count+1
              # Source Cone
               elif (middle < end1) and (middle < end2):
                       source_cone[imlt][iL]=source_cone[imlt][iL]+1
                       plot_title='Source Cone'
                       count=count+1
              # Loss Cone
               elif (middle > end1) and (middle > end2):
                       loss_cone[imlt][iL]=loss_cone[imlt][iL]+1
                       plot_title='Loss Cone'
                       count=count+1
              # One sided cone
               elif (end1 < middle) and (middle < end2) or (end2 < middle) and (middle < end1):
                      one_cone[imlt][iL]=one_cone[imlt][iL]+1
                      plot_title='One Sided Source'
                      count=count+1
              # Now for the cases where the endpoints are compromised, put these into uncategorized
              # Source cone
               elif (middle < sec1) and (middle < sec2):
                       source_cone[imlt][iL]=source_cone[imlt][iL]+1
                       plot_title='Source Cone'
                       count=count+1
              # Loss Cone
               elif (middle > sec1) and (middle > sec2):
                       loss_cone[imlt][iL]=loss_cone[imlt][iL]+1
                       plot_title='Loss Cone'
                       count=count+1
              # One sided cone, leave for compromised end points for now
               elif (sec1 < middle) and (middle < sec2) or (sec2 < middle) and (middle < sec1):
                      one_cone[imlt][iL]=one_cone[imlt][iL]+1
                      plot_title='One Sided Source'
                      count=count+1  
              # doesn't match anything
               else:
                      uncat[imlt][iL]=uncat[imlt][iL]+1
                      count=count+1
                      plot_title='Uncategorized'

    
              # make a quick plot of this and save it
             # time_label=datetime.datetime.strftime(epoch[itime+window/2], '%d %b %Y %H:%M:%S')
             # tl2=datetime.datetime.strftime(epoch[itime+window/2],'%Y%m%d')
             # fig=plt.figure()
             # ax=fig.add_subplot(111)
             # fig.subplots_adjust(bottom=0.15, left=0.15, right=0.9)
             # ax.set_title(plot_title + ' at '+ time_label, fontsize=18, fontweight='bold')
             # ax.set_ylabel('Normalized Flux', fontsize=22, fontweight='bold')
             # ax.set_xlabel(r'PA Bins', fontsize=22, fontweight='bold')
             # ax.set_yscale('log')
             # ax.set_xlim(0, 180)
             # font = {'family' : 'normal',
             #         'weight' : 'bold',
             #         'size'   : 22}
             # plt.rc('font', **font)
             # ax.scatter(PA_labels, nflux, c='b', s=100,marker='x', lw=3)
             # ax.scatter(np.array(PA_labels)[Low_indices], nflux[Low_indices], c='DarkOrange', s=50, marker='o') 
             ## ax.plot(np.linspace(-1, 1, 181),fit, c='r', lw=3)
             # ax.plot(mid_points , [end1, sec1, middle, sec2, end2], c='r', lw=3)
             # ax.plot(np.linspace(0, 180, 181),np.ones(181)*peak_value, lw=3, ls='--', c='g')
             # ax.plot(np.linspace(0, 180, 181),np.ones(181)*min_value, lw=3, ls='--', c='g')
             # ax.annotate( 'L = ' + str(iL*0.25) + '\n MLT = '+ str(imlt*0.5),xy=(.95,0.95), xycoords='axes fraction',  horizontalalignment='right', verticalalignment='top', fontsize=20)
             # os.chdir('All_PA_Combine_Plots')
             # subdir_name=plot_title
             # if not os.path.exists(subdir_name):
             #     os.umask(0) # unmask if necessary
             #     os.makedirs(subdir_name, 0777) 
             # os.chdir(subdir_name)#
             # plt.savefig('combined_'+str(tl2)+'_'+plot_title+'_'+str(count)+'_n='+str(ien)+'.pdf')
             # plt.close()
             # os.chdir('..')
             # os.chdir('..')
               break
# add a day
        dt1=dt1+datetime.timedelta(days=1)
        cEn.append(np.nanmean(np.array(temp_en)))
  # F_lkp gets pretty big by the end
  # now pickle it for the average energy channel
     avg_en=int(np.nanmean(np.array(cEn))*100)/100.0
     print cur_date+ '_energy='+str(avg_en)
    # #
     os.chdir('PA_Combine_Bins/Window='+str(window))
     cats=[inverse_butterfly, source_cone, one_cone, butterfly, isotropic, loss_cone,uncat, bad_counts]
     for icat in range(len(cats)):
        with open(cur_date+'_'+species[ispe]+'_binned_fluxes_combine_'+ name[idir]+'_'+str(ien)+'_'+str(avg_en)+'_cat='+str(cat_labels[icat])+'.pickle', 'wb') as f:
            pickle.dump(cats[icat], f)
     with open(cur_date+'_'+species[ispe]+'_binned_fluxes_combine_'+ name[idir]+'_'+str(ien)+'_'+str(avg_en)+'_missing_fluxes.pickle', 'wb') as f2:
            pickle.dump(missing_fluxes, f2) # should dump as a 3-D array 11 x Lbins X MLTbins
     os.chdir('..')
     os.chdir('..')


        



