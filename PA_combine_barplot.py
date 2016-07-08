# PA_combine_barplot.py
#
# plot the totals of all the sections
#
# LKS, March 2015
#
# imports 
import numpy as np
import os
import pickle
import datetime
import matplotlib as mpl
from matplotlib import pyplot as plt
import itertools as itert
from matplotlib.colors import LogNorm, BoundaryNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import glob
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
from dateutil.relativedelta import relativedelta 
#  
class HandlerSquare(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = xdescent + 0.5 * (width - height), ydescent
        p = mpatches.Rectangle(xy=center, width=height,
                               height=height, angle=0.0)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]   
#
# globals
date1='20130201'
date2='20150401'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
color_array=['r', 'b', 'g']
def diff_month(d1, d2):
    return (d1.year - d2.year)*12 + d1.month - d2.month

window=10
total_month=diff_month(dt2, dt0)
species=['H']
directory=['rbspa', 'rbspb']
cats=[ 'butterfly', 'inverse_butterfly', 'isotropic', 'loss_cone', 'one_cone', 'source_cone', 'uncategorized', 'bad_counts','spacecraft_charging']
mltbins=48
mltbinsR=24
Lbins=11
LBinsR=5
#
# get the files
for ispe in range(len(species)):
    for ien in range(0,16):
        uncat=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        inv_but=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        s_cone=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        o_cone=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        but=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        iso=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        l_cone=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        sc_charging=[ [ 0 for x in xrange(Lbins) ] for x in xrange(mltbins)]
        data_arr=[[ np.array(0) for x in xrange(mltbinsR) ] for x in xrange(len(cats))]
        data_total=[ 0  for x in xrange(mltbinsR)]
        missing_fluxes=[[[ 0 for x in xrange(Lbins)] for x in xrange(mltbins)] for x in xrange(11)]

        dt1=dt0    
        for imonth in range(total_month):

         cur_date=str(dt1.month)+'_'+str(dt1.year)
         dt1=dt1+relativedelta(months=1)
         os.chdir('/Users/loisks/Documents/ResearchProjects/SummedPitchAngles/Window10')
         for icat in range(len(cats)):
          get_files_a=glob.glob(cur_date+'_'+species[ispe]+'_sc_binned_fluxes_combine_rbspa_'+str(ien)+'_'+'*')
          get_files_b=glob.glob(cur_date+'_'+species[ispe]+'_sc_binned_fluxes_combine_rbspb_'+str(ien)+'_'+'*')
          for iline in range(len(get_files_a)):
              if cats[icat] in get_files_a[iline]:
                  index=iline
                  with open(get_files_a[index], 'rb') as f:
                      data1=pickle.load(f)
                  with open(get_files_b[index], 'rb') as f2:
                      data2=pickle.load(f2)
# WHEN I GET BACK, check to make sure data is loading and we are getting real values
                  for imlt in range(mltbinsR):
                      for iL in range(2,Lbins-5): # L between 2-3
                          data_arr[icat][imlt]+= data1[imlt*2][iL] + data2[imlt*2][iL] + data1[(imlt*2)+1][iL] + data2[imlt*2+1][iL] 
              else:
                  continue
         os.chdir('..')
         os.chdir('..')
        for icat in range(len(cats)):
            for imlt in range(mltbinsR):
                    data_total[imlt]+=data_arr[icat][imlt]
 #
#+'_'
        # for imlt in range(mltbins):
        #      for iL in range(Lbins):
        #          total[imlt][iL]= len(data_total[imlt][iL])
        print "total number of PADs " + str(np.sum(data_total))
        print "going into plot"
        for icat in range(len(cats)):
              temp_data=data_arr[icat] # now its MLT 
              temp_total=data_total
              plot_title=cats[icat]
              fig=plt.figure()
              L_cur=" 2 to 3"
              Lcur="2_to_3"
              fig.subplots_adjust(bottom=0.15, left=0.16, right=0.9)
              font = {'family' : 'normal',
                      'weight' : 'bold',
                      'size'   : 22}
              plt.rc('font', **font)
              # combine the bad counts and uncategorized
              ax=fig.add_subplot(111)
              #ax.set_title(plot_title + ' at L=' + str(L_cur))
              ax.set_ylabel('# of PADs', fontweight='bold')
              ax.set_xlabel('MLT', fontweight='bold')
              ax.set_xlim(0, 24)
              ax.set_ylim(0,3000)

              try:
                ax.bar(range(mltbinsR), temp_total, align='center', color='gray', alpha=0.5)
                ax.bar(range(mltbinsR),temp_data, align='center', color=color_array[ispe])
          #    except(ValueError):
              except(IOError):
                  print icat
                  continue
              label_array=[]; count=0
              for imlt in range(6):
                  label_array.append(str(imlt*4))
                  label_array.append(' ')
                  label_array.append(' ')
                  label_array.append(' ')
              ax.set_xticks(range(mltbinsR))
              ax.set_xticklabels(label_array)
              os.chdir('/Users/loisks/Documents/ResearchProjects/SummedPitchAngles/SC_BarPlots/Window10')
              plt.savefig('Combine_sc_'+species[ispe]+'energy='+str(ien)+'L='+str(Lcur)+'_'+plot_title+'n='+str(ien)+'.pdf')
              plt.close()
              os.chdir('..')
              os.chdir('..')
#
#
# STACKED BAR
        fig=plt.figure()
        L_cur=" 2 to 3"
        Lcur="2_to_3"
        fig.subplots_adjust(bottom=0.15, left=0.16, right=0.79)
        font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}
        plt.rc('font', **font)
        # combine the bad counts and uncategorized
        plot_title='stacked_all'
        ax=fig.add_subplot(111)
        #ax.set_title(plot_title + ' at L=' + str(L_cur))
        ax.set_ylabel('# of PADs', fontweight='bold')
        ax.set_xlabel('MLT', fontweight='bold')
        ax.set_xlim(0, 24)
        ax.set_ylim(0,3000)
        loss=np.array(data_arr[3]); one=np.array(data_arr[4])
        source=np.array(data_arr[5]); bad=np.array(data_arr[7])
        
        r1=ax.bar(range(mltbinsR), data_total, color='gray', alpha=0.5)
        r2=ax.bar(range(mltbinsR),loss, color='DarkOrange') # loss cone
        r3=ax.bar(range(mltbinsR),one, color='Red', bottom=loss) # one cone
        r4=ax.bar(range(mltbinsR),source, color='Yellow', bottom=loss+one) # source cone
        r5=ax.bar(range(mltbinsR),bad, color='k', bottom=loss+source+one) #  bad
        params = {'legend.fontsize': 15}
        plt.rcParams.update(params)
        plt.legend([r1[0], r2[0], r3[0], r4[0], r5[0]], ['Other','Loss', 'One Side', 'Source', 'Uncat.'], bbox_to_anchor=[1.35, 0.75] 
                )
        
        
        label_array=[]; count=0
        for imlt in range(12):
            label_array.append(str(imlt*4))
            label_array.append(' ')
            label_array.append(' ')
            label_array.append(' ')
        ax.set_xticks(range(mltbinsR))
        ax.set_xticklabels(label_array)

 
        os.chdir('/Users/loisks/Documents/ResearchProjects/SummedPitchAngles/SC_BarPlots/Window10')#
        plt.savefig('updated_lo_sc_Combine_'+species[ispe]+'energy='+str(ien)+'L='+str(Lcur)+'_'+plot_title+'n='+str(ien)+'.pdf')
        plt.close()
        os.chdir('..')
        os.chdir('..')
#
#
# STACKED BAR PERCENTAGE
#
#
        fig=plt.figure()
        L_cur=" 2 to 3"
        Lcur="2_to_3"
        fig.subplots_adjust(bottom=0.15, left=0.16, right=0.79)
        font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}
        plt.rc('font', **font)
        # combine the bad counts and uncategorized
        plot_title='percentage'
        ax=fig.add_subplot(111)
        #ax.set_title(plot_title + ' at L=' + str(L_cur))
        ax.set_ylabel('% of Total PADs', fontweight='bold')
        ax.set_xlabel('MLT', fontweight='bold')
        ax.set_xlim(0, 24)
        ax.set_ylim(0,100)
        data_total=np.array(data_total)
        CountsSummed=np.sum(data_total)
        UncatPercent=np.sum(data_arr[7])/(1.0*CountsSummed)
        print UncatPercent
        print 'n = ' + str(ien)
        lossP=(1.0*np.array(data_arr[3])/data_total)*100
        oneP=(1.0*np.array(data_arr[4])/data_total)*100
        sourceP=(1.0*np.array(data_arr[5])/data_total)*100
        badP=(1.0*np.array(data_arr[7])/data_total)*100
        #        
        r1=ax.bar(range(mltbinsR), (data_total/data_total)*100, color='gray', alpha=0.5)
        r2=ax.bar(range(mltbinsR),lossP,  color='DarkOrange') # loss cone
        r3=ax.bar(range(mltbinsR),oneP, color='Red', bottom=lossP) # one cone
        r4=ax.bar(range(mltbinsR),sourceP, color='Yellow', bottom=lossP+oneP) # source cone
        r5=ax.bar(range(mltbinsR),badP, color='Black', bottom=lossP+sourceP+oneP) #  bad
        params = {'legend.fontsize': 15}
        plt.rcParams.update(params)
        plt.legend([r1[0], r2[0], r3[0], r4[0], r5[0]], ['Other','Loss', 'One Side', 'Source', 'Uncat.'], bbox_to_anchor=[1.35, 0.75] 
                )
        label_array=[]; count=0
        for imlt in range(6):
            label_array.append(str(imlt*4))
            label_array.append(' ')
            label_array.append(' ')
            label_array.append(' ')
        ax.set_xticks(range(mltbinsR))
        ax.set_xticklabels(label_array)
        os.chdir('SC_BarPlots/Window'+str(window))#
        plt.savefig('updated_sc_lo_Combine_'+species[ispe]+'energy='+str(ien)+'L='+str(Lcur)+'_'+plot_title+'n='+str(ien)+'.pdf')
        plt.close()
        os.chdir('..')
        os.chdir('..')
                    
