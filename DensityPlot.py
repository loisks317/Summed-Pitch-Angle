# DensityPlot.py
#
# combines satellite A and B data 
# plots the spacecraft potential corrected density
# uses the pickling function for hdf5_data_save
#
# LKS, November 2015 -- First Snow of the Year --
#
#

# imports 
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
os.chdir('/Users/loisks/Desktop/Functions/')
import pickling
import plots
os.chdir('/Users/loisks/Documents/ResearchProjects/SummedPitchAngles/')


nLbins=11
nmlt_bins=48
mData=[ [[] for i in range(nLbins)] for j in range(nmlt_bins)]
#
# load in files and plot
os.chdir('CorrectedDensity')
data=pickling.hdf5_data_open('combined_plasma_density_0_16_species=P.h5', nLbins, nmlt_bins)
os.chdir('..')
for i in range(nmlt_bins):
    for j in range(nLbins):
        mData[i][j]=np.nanmedian(data[i][j])
        #if np.isnan(mData[i][j])==True:
        #    mData[i][j]=0
plots.fancy_plot(mData, nmlt_bins, nLbins, -1, 2, 'Density [cm$^{-3}$]', 'H$^{+}$', 'kp3', 'HDensity', 'CorrectedDensity', 'log', '')
