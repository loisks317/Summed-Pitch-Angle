# PitchAngleClassification.py
#
# This code sums pitch angle distributions by PA Bin in HOPE
# We specifically look at 1-10 eV here
# Classify the shape of the distribution using pre-set criteria
# 2015-12-7, now includes spacecraft potential
# 
#
#
# LKS March 2015, edited June 2015 & November 2015
# For publication in GRL, Sarno-Smith 2015
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
import matplotlib.pyplot as plt
import h5py
from scipy import interp
from dateutil.relativedelta import relativedelta
#
# INPUT PARAMETERS
#
# start date
date1='20130201'
#
# end date
#date2='20130401'
date2='20150401'
#
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
    # gives the number of months in time range for easy file output
    return (d1.year - d2.year)*12 + d1.month - d2.month
total_month=diff_month(dt2, dt0)
#
# number of LShell Bins and MLT bins
# MLT range is 0 to 48, with 0.5 increments
# L-Shell range is 1.5 to 4, with 0.25 increments
Lbins=11
mltbins=48
Lmax=4
Lmin=1.5

#
# set limitations
kpBar=3 # threshold for Kp
countsBar=10 # how many counts required in a 'good' PA bin
binBar=5 # number of bins with enough counts required for PAD
HopeEnergies=[.99,1.2,1.34,1.55,1.83,2.18,2.53,2.95,3.38,3.94,4.64,5.35,6.26,7.32,8.51,9.92]
# energy channels
n1=0 # 0.99 eV
n2=16 # 9.9 eV
PAs=11 # 11 pitch angles on HOPE

#
# bins to sum
window=10
PA_labels=[4.5, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0, 175.5]
mid_points=np.array(  [13.5, 45, 90, 135, 166.5])

#
# Satellite Directories and Species types
dir=['HOPE_PA_A', 'HOPE_PA_B']
dirE=['EFW_L3_A', 'EFW_L3_B']
name=['rbspa', 'rbspb']
species=['H']#, 'He', 'O']
name_species=[ 'FPDU']#, 'FHEDU','FODU'] # gives the unidirectional flux with PA

os.chdir('/Users/loisks/Desktop/PA_check/PA_SC_Combine_Bins/')
if not os.path.exists('Window'+str(window)):
        os.umask(0) # unmask if necessary
        os.makedirs('Window'+str(window), 0777) 
os.chdir('..')
#
# Loop over Species
for ispe in range(len(name_species)):
#
# Loop over Energy    
 for ien in range(n1,n2):
#
# Loop over directory
  for idir in range(len(dir)):
   dt1=dt0 # start over at the beginning date
#
# Loop over months
   for imonths in range(total_month):

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
    allBad=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    countsHisto={}
    # set up counts histogram
    for iC in range(PAs+1):
        countsHisto[iC]=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    uncat=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    spacecraft_charging=[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)]
    missing_fluxes=[[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)] for x in xrange(len(PA_labels))]
# Categories
    cats=[inverse_butterfly, source_cone, one_cone, butterfly, isotropic, loss_cone, uncat, bad_counts, spacecraft_charging]
    cat_labels=[ 'inverse_butterfly', 'source_cone', 'one_cone', 'butterfly', 'isotropic', 'loss_cone', 'uncategorized','bad_counts', 'spacecraft_charging']
#
# See if we already sorted for it
# CHANGE TO SC LATER
    #os.chdir('/Users/loisks/Desktop/PA_check/PA_SC_Combine_Bins/Window10/')
    #check=glob.glob(cur_date+'_'+species[ispe]+'_sc_binned_fluxes_combine_'+ name[idir]+'_'+str(ien)+'*')
    #try:
    #    check[0] # see if anything is there ,this will result in an index error if nothing is there
    #    print cur_date+'_'+species[ispe]+'_'+name[idir]+'_'+str(ien)
    #    dt1=dt2
# if nothing is there, proceed to loop through days in month
    #except(IndexError):
    while dt1 < dt2:
       while True:
            try:
                date1=datetime.datetime.strftime(dt1, '%Y%m%d')
                # put link where you store files here
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/'+dir[idir])
                hope_file=name[idir]+'_rel03_ect-hope-PA-L3_'+date1+'*.cdf'
                ghope=glob.glob(hope_file)
                pyf=pycdf.CDF(ghope[0])
                #
                # now get s/c potential data
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/'+dirE[idir])
                EFWfile='*'+date1+'*.cdf'
                f=glob.glob(EFWfile)
                pyfE=pycdf.CDF(f[0])
                break
            except:
                print dt1
                dt1=dt1+datetime.timedelta(days=1)
                pass

        # data quantities from HOPE
       try:
        L=pyf['L_Ion'][...] # L Shell
        energy=pyf['HOPE_ENERGY_Ion'][...] # Hope Energy
        MLT=pyf['MLT_Ion'][...] # MLT
        FLUX=pyf[name_species[ispe]][...] # Differential number Flux
        COUNTS=pyf['Counts_P'][...] # time x 11 x 72
        epoch=pyf['Epoch_Ion'][...]    
        PA_LABEL=pyf['PITCH_ANGLE'][...] # right PA label
        EFWepoch=pyfE['epoch'][...]
        EFWSC=-1*np.array(pyfE['Vavg'][...])
        os.chdir('/Users/loisks/Desktop/PA_check')
#
# Use spacepy to get the kyoto Kp data
# spacepy is available at http://spacepy.lanl.gov/
# I highly recommend it for loading cdf files in python
# and fetching things like Kp and Dst
# thanks SpacePy team!
        kyoto=spk.fetch('kp',dt1, dt1)
        day=datetime.datetime.strftime(dt1,'%d') # need to get day of month
        kp_arr=kyoto['kp'][(int(day)-1)*8: int(day)*8]
        # get the kp every 3 hours
        hours=np.array([0, 3, 6, 9, 12, 15, 18, 21])
        #
        temp_en=[]
        count=0
        for itime in range(len(epoch)/window): 
            while True: # everything ends in a break
               # gives the L indexes
               iLArr=np.array(L[itime*window:(window*itime+window)])/0.25
               iL=int(np.round(np.nanmean(iLArr))-(Lmin/0.25))
               #if (iL*4+Lmin <= 2*4) or (iL*4+Lmin>=3*4): # if this is greater than zero # only consider 2 < L < 3
               if (iL >  Lbins-1):
                   break # this means that we're beyond our area of interest
               iMltArr=np.array(MLT[itime*window:(window*itime+window)])/0.5
               imlt=int(np.round(np.nanmean(iMltArr)))
               if imlt ==48: # as it should be
                   imlt=0
               cengArr=np.swapaxes(np.array(energy[itime*window:(window*itime+window)]), 1, 0)[n1:n2]# 10 x 16 arrayarray
               cengArr=np.nanmean(cengArr, axis=1) # now just 16 elements
               if len(np.where(cengArr > 20)[0]) > 0:
                   # perigee mode for HOPE, not good for our 1-10 eV study 
                   break
               # Kp!
               middleTime=epoch[window*itime+(window/2)]
               # do the same approach to get s/c potential
               hourIndex=np.where(hours <= middleTime.hour)[0][-1]
               if kp_arr[hourIndex] > kpBar:
                   break # high activity
               try:
                   phiIbeg=np.where(EFWepoch >= epoch[window*itime])[0][0]
               except(IndexError):
                   phiEbeg=0
               try:
                   phiIend=np.where(EFWepoch >= epoch[window*itime+window])[0][0]
               except(IndexError):
                   phiIend=-1
               try:
                  cEnArr=np.array(cengArr)+np.abs(np.nanmedian(EFWSC[phiIbeg:phiIend])) #uniformly add the potential
               except:
                   break
               cfluxArr=np.array(np.swapaxes(FLUX[(window*itime):(itime*window+window)], 2, 0)) # ready for ien
               cfluxArr[cfluxArr > 1e12]=np.nan # bad data. 
               cfluxArr[cfluxArr < 0]=np.nan # also bad data.
               cCountsArr=np.array(np.swapaxes(COUNTS[(window*itime):(itime*window+window)], 2, 0))
               
               ENAVG=HopeEnergies[ien]
               try:
                      EnArr=np.array(cEnArr)[~np.isnan(cEnArr)]
                      below=np.max(EnArr[EnArr<=HopeEnergies[ien]])
                      BelowI=np.where(EnArr==below)[0]
                      above=np.min(EnArr[EnArr>=HopeEnergies[ien]])
                      AboveI=np.where(EnArr==above)[0]
                      f=(HopeEnergies[ien]- below)/(1.0*((HopeEnergies[ien]-below)+(above-HopeEnergies[ien])))
                      intpFlux=(np.array(cfluxArr[AboveI])**f)*(np.array(cfluxArr[BelowI])**(1-f))[0]
                      intpCounts=(np.array(cCountsArr[AboveI])**f)*(np.array(cCountsArr[BelowI])**(1-f))[0]
                      # check to make sure both of these are 2D
                      # this gives a result for a specific hope energy channel so we don't need to change binning :D
                      intpFlux=intpFlux[0]
                      intpCounts=intpCounts[0]
                      intpFlux[intpFlux==0]=np.nan
                      print 'nonzero flux array'

               except(ValueError): # when there is no flux b/c of s/c charging
                       intpFlux=np.array([[np.nan]*11]*window)
                       intpCounts=np.array([[np.nan]*11]*window)
                
               
               #temp_en.append(cengArrMean)
               # of the flux based on the counts
               summedCounts=np.nansum(intpCounts, axis=1)
               notValidCounts=np.where(summedCounts < countsBar)[0] # define a hard limit
               # 
               # do a weighted average of the fluxes
               weightedMultiplication = np.array(intpFlux)*np.array(intpCounts)
               # now sum along the time axis
               summedWeights=np.nansum(weightedMultiplication, axis=1)
               weightedFlux=summedWeights/summedCounts
               weightedFlux[notValidCounts]=np.nan
               #
               # rarely will have a full distribution now
               # due to counting issues
               # INSERT HISTOGRAM COUNTER HERE
               Low_indices=np.where(np.isnan(weightedFlux)==True)[0]

               #

               for iCount in range(PAs+1):
                   if len(Low_indices) == iCount: # number of bad indicies
                       countsHisto[iCount][imlt][iL]+=1
                       break # get out of this loop
                   
               # requirement on number of bins
               if len(Low_indices) > binBar:
                  bad_counts[imlt][iL]+=1
                  if (len(Low_indices) == 11):
                      allBad[imlt][iL]+=1
                    #  print "sc charging all zeros at " + str(epoch[itime+window/2])
                  break # PAD is mostly zeros, just get rid of it

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
               nflux=(weightedFlux/nfactor) # now it's normalized
               zeros=len(np.where(nflux==0)[0])
               nans=len(np.where(np.isnan(nflux)==True)[0])
               if (zeros + nans) > binBar: # bad data
                  # print "spacecrafting charging event"
                   spacecraft_charging[imlt][iL]+=1
                   print "sc charging nans and zeros at " + str(epoch[itime+window/2])
                   break

               # order the normalized flux without the nans
               oflux = np.sort(nflux[~np.isnan(nflux)])
               min_value=oflux[1] # exclude outlier 
               peak_value=oflux[-2] # exclude outlier

               # check for isotropy
               isoCheck=peak_value/min_value
               #
               # get the means in each section
               end1=np.nanmean(nflux[0:2])
               sec1=np.nanmean(nflux[2:4])
               middle=np.nanmean(nflux[4:7])
               sec2=np.nanmean(nflux[7:9])
               end2=np.nanmean(nflux[9:11])
               
               # 
               # ALGORITHM PART
               #
               
               #
               # Isotropic
               if isoCheck < 1.2: # if the maximum is less than 1.2x the minimum
                    isotropic[imlt][iL]=isotropic[imlt][iL]+1
                    plot_title='Isotropic'
                    count=count+1
                    print weightedFlux
                    print 'iso sort'
               #
               # Butterfly
               elif (end1 < sec1) and (end2 < sec2) and (sec1 > middle) and (sec2 > middle):
                butterfly[imlt][iL]=butterfly[imlt][iL]+1
                plot_title='Butterfly'
                count= count+1
                print 'butterfly sort'

              #  
              # Inverse Butterfly
               elif (end1 > sec1) and (end2 > sec2) and (sec1 < middle) and (sec2 < middle):
                  inverse_butterfly[imlt][iL]=inverse_butterfly[imlt][iL]+1
                  plot_title='Inverse Butterfly'
                  count= count+1
                  print 'inv. butterfly sort'

              #
              # Source Cone
               elif (middle <= end1) and (middle <= end2):
                       source_cone[imlt][iL]=source_cone[imlt][iL]+1
                       plot_title='Source Cone'
                       count=count+1
                       print  'source sort'
              #
              # Loss Cone
               elif (middle > end1) and (middle > end2):
                       loss_cone[imlt][iL]=loss_cone[imlt][iL]+1
                       plot_title='Loss Cone'
                       count=count+1
                       print 'loss sort'

              #         
              # One sided cone
               elif (end1 <= middle) and (middle <= end2) or (end2 <= middle) and (middle <= end1):
                      one_cone[imlt][iL]=one_cone[imlt][iL]+1
                      plot_title='One Sided Source'
                      count=count+1
                      print 'one side sort'

              #
              # Now for the cases where the endpoints are compromised, put these into uncategorized
              #
              
              # Compromised Source cone
               elif (middle < sec1) and (middle < sec2):
                       source_cone[imlt][iL]=source_cone[imlt][iL]+1
                       plot_title='Source Cone'
                       count=count+1
                       print 'source'
              #
              # Compromised Loss Cone
               elif (middle > sec1) and (middle > sec2):
                       loss_cone[imlt][iL]=loss_cone[imlt][iL]+1
                       plot_title='Loss Cone'
                       count=count+1
                       
                       print 'loss cone'
              #
              # One sided cone, leave for compromised end points for now
               elif (sec1 < middle) and (middle < sec2) or (sec2 < middle) and (middle < sec1):
                      one_cone[imlt][iL]=one_cone[imlt][iL]+1
                      plot_title='One Sided Source'
                      count=count+1
                      print 'one side'
              # doesn't match anything
               else:
                      uncat[imlt][iL]=uncat[imlt][iL]+1
                      count=count+1
                      plot_title='Uncategorized'
                      print sec1, middle, sec2, end1, end2
                      print weightedFlux

             # this makes the lineplot of Figure 2
             # make a quick plot of this and save it
               time_label=datetime.datetime.strftime(epoch[itime+window/2], '%d %b %Y %H:%M:%S')
               tl2=datetime.datetime.strftime(epoch[itime+window/2],'%Y%m%d')
               fig=plt.figure()
               ax=fig.add_subplot(111)
               fig.subplots_adjust(bottom=0.15, left=0.15, right=0.9)
               ax.set_title(plot_title + ' at '+ time_label, fontsize=18, fontweight='bold')
               ax.set_ylabel('Normalized Flux', fontsize=22, fontweight='bold')
               ax.set_xlabel(r'PA Bins', fontsize=22, fontweight='bold')
               ax.set_yscale('log')
               ax.set_xlim(0, 180)
               font = {'family' : 'normal',
                       'weight' : 'bold',
                       'size'   : 22}
               plt.rc('font', **font)
               ax.scatter(PA_labels, nflux, c='b', s=100,marker='x', lw=3)
               ax.scatter(np.array(PA_labels)[Low_indices], nflux[Low_indices], c='DarkOrange', s=50, marker='o') 
            #   ax.plot(np.linspace(-1, 1, 181),fit, c='r', lw=3)
               ax.plot(mid_points , [end1, sec1, middle, sec2, end2], c='r', lw=3)
               ax.plot(np.linspace(0, 180, 181),np.ones(181)*peak_value, lw=3, ls='--', c='g')
               ax.plot(np.linspace(0, 180, 181),np.ones(181)*min_value, lw=3, ls='--', c='g')
               ax.annotate( 'L = ' + str(iL*0.25) + '\n MLT = '+ str(imlt*0.5),xy=(.95,0.95), xycoords='axes fraction',  horizontalalignment='right', verticalalignment='top', fontsize=20)
               os.chdir('All_PA_Combine_Plots')
               subdir_name=plot_title
               if not os.path.exists(subdir_name):
                   os.umask(0) # unmask if necessary
                   os.makedirs(subdir_name, 0777) 
               os.chdir(subdir_name)#
               plt.savefig('combined_'+str(tl2)+'_'+plot_title+'_'+str(count)+'_n='+str(ien)+'.pdf')
               plt.close()
               os.chdir('..')
               os.chdir('..')
               break
# add a day
        dt1=dt1+datetime.timedelta(days=1)
        #cEn.append(np.nanmean(np.array(temp_en)))
       except(OSError):
           #try again
           pass
  # F_lkp gets pretty big by the end
  # now pickle it for the average energy channel
#     avg_en=HopeEnergies[ien]
#     print cur_date+ '_energy='+str(avg_en)
#    # #
#     os.chdir('PA_SC_Combine_Bins/Window'+str(window))
#     cats=[inverse_butterfly, source_cone, one_cone, butterfly, isotropic, loss_cone,uncat, bad_counts]
#     print cats[5]
#     for icat in range(len(cats)):
#        with open(cur_date+'_'+species[ispe]+'_sc_binned_fluxes_combine_'+ name[idir]+'_'+str(ien)+'_'+str(avg_en)+'_cat='+str(cat_labels[icat])+'.pickle', 'wb') as f:
#            pickle.dump(cats[icat], f)
#     with open(cur_date+'_'+species[ispe]+'_sc_binned_fluxes_combine_'+ name[idir]+'_'+str(ien)+'_'+str(avg_en)+'_missing_fluxes.pickle', 'wb') as f2:
#            pickle.dump(missing_fluxes, f2) # should dump as a 3-D array 11 x Lbins X MLTbins
#     for iC in range(12):
#         pickle.dump( np.array(countsHisto[iC]),open(cur_date+'_sc_counts='+str(iC)+'_'+species[ispe]+'_'+str(ien)+'_'+name[idir]+'_'+str(avg_en)+'_countsHisto.p', 'wb'))
#     os.chdir('..')
#     os.chdir('..')


        



