#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 13:28:38 2021

@author: laurahenry
"""
import math 
import scipy 
import glob # for listing files in a directory
import os # for sep, ...
import os.path
import numpy as np
import pylab
import matplotlib.pyplot as plt
from pylab import *
pylab.ion() # switch on interactive figure mode?
from scipy import interpolate

#Extension of the last file to process 
filenumber = 47

#Filename (everything before 'chxx_00xx.x_y') and directory where your data is stored
filename='MgO_Au_'
filepath = ''
# Y gap Offset for stacked data in figure(1) 
offset = 1   

# Calibration parameters (this needs to be updated by the user for the specific experiment using FalxconXCalib_....txt)

Calib_E1 = [0.049779162, 0.048955129, 0.049851929, 0.051011480, 0.050300672, 0.049917061, 0.049435997]
Calib_E2 = [0.00000014, 0.00000013, 0.00000013, 0.00000013, 0.00000014, 0.00000006, 0.00000014]
theta=[8.103, 8.045, 8.02, 8.006, 8.029, 8.065, 8.125]


#Begin of code for combining all 7 channel into one 2D integrated spectrum. 
channels=np.arange(0,7,1)
DATA=np.zeros((len(channels),2047,filenumber))
real_filenumber = 0
x=np.zeros([7,2047])
xd = np.arange(1, 3, 2/2047)  #used for d-spacing presentation (lower, higher, (higher-lower)/2047)
xE = np.arange(30, 70, 40/2047)   #used for Energy presentation  (lower, higher, (higher-lower)/2047)
xQ = np.arange(2, 5.5, 3.5/2047)   #used for Energy presentation  (lower, higher, (higher-lower)/2047)


def Convert_E_channel(E1,i):
    Chan=int(-(-Calib_E1[i] + math.sqrt(Calib_E1[i]**2 - 4*Calib_E2[i]*E1))/(2*Calib_E2[i]))
    return Chan
def Convert_channel_E(chan,i):
    E=float(chan*Calib_E1[i]+chan**2*Calib_E2[i])
    return E
    return Chan
def Convert_E_dsp(E,i):
    dsp=299792458.*6.62607E-34/2/sin(theta[i]*pi/180./2.)/E*0.001/1.60217E-19*10000000000.
    return dsp   
def Convert_dsp_Q(dsp):
    Q=2*pi/dsp
    return Q  


for i in range(len(x)):
    for j in range(2047):
        x[i,j]=j


for i, j in enumerate(channels):
    for jj in range(2047):
        x[i,jj]=Convert_channel_E(x[i,jj], i)  #Default is energy if below is commented (+++ change line 97 for xE)
        x[i,jj]=Convert_E_dsp(x[i,jj], i) #uncomment  with above for d-spacing axis (+++ change line 97 for xd)
        # x[i,jj]=Convert_dsp_Q(x[i,jj]) #uncomment with above for Q range (+++ change line 97 for xQ)


for i,j in enumerate(channels):
    for ii in range(filenumber):
        if os.path.isfile(filepath+ filename+"ch0%s_000%s.x_y" %(j, ii+1)) or os.path.isfile(filepath+ filename+"ch0%s_00%s.x_y" %(j, ii+1)):
            if ii < 9:
                data = np.loadtxt(open(filepath+filename+"ch0%s_000%s.x_y" %(j, ii+1), 'rb') ,delimiter =' ',skiprows=1)
                print(filename+"ch0%s_000%s.x_y" %(j, ii))
            else:
                data = np.loadtxt(open(filepath+filename+"ch0%s_00%s.x_y" %(j, ii+1), 'rb') ,delimiter =' ',skiprows=1)
                print(filename+"ch0%s_00%s.x_y" %(j, ii+1))
            DATA[i,:,ii] = data[:,1]
            real_filenumber +=1
        else:
            print('There is no such file')
real_filenumber//=len(channels)
DATA2 = np.asarray(DATA)
#plot(DATA[0,:,:])
DATA = np.zeros((len(channels),2047,real_filenumber))
DATA_combined = np.zeros([2047,real_filenumber])
DATA_combined = np.asarray(DATA_combined)

file=0
for ii in range(filenumber):
    if DATA2[:,:,ii].max()>0:
        DATA [:,:,file] += DATA2[:,:,ii]   #Raw intensities
        DATA [:,:,file] += DATA2[:,:,ii]/DATA2[:,:,ii].max()  #Normalized Intensity
        file+=1
    else:
        ii+=1
ii=0
        
        
for i,j in enumerate(channels):
    for ii in range(real_filenumber):
        interpol1=interpolate.interp1d(x[i,:],DATA[i,:,ii], kind='linear',bounds_error=False, fill_value="extrapolate")
        DATA[i,:,ii] = interpol1(xd)   ### options are xE, xQ, xd
            
colors = plt.cm.viridis(np.linspace(0,1,filenumber))
for ii in range(real_filenumber):
    DATA_combined[:,ii]=DATA[:,:,ii].sum(0)
    figure(2)
    
    """Plot as a function of d-spacing"""
    plot(xd, DATA_combined[:,ii]+ii*offset,color = colors[ii])
    xlabel('d-spacing (Angstrom)') 
    index = 'd-spacing'

    """Plot as a function of E"""
    # plot(xE, DATA_combined[:,ii]+ii*offset,'k')
    # xlabel('Energy (keV)') 
    #index = 'Energy (keV)'



    """Plot as a function of Q"""
    #plot(xQ, DATA_combined[:,ii]+ii*offset,'k')
    #xlabel('Q (Angstrom -1)') 
    #index = 'Q (Angstrom -1)'

    
    ylabel('Intensity (a.u.)')
    
DATA_combined = np.asarray(DATA_combined)

np.savetxt(filepath + str(filename)+'_datas.txt', DATA_combined, delimiter=' ')
np.savetxt(filepath + str(filename)+'d_spacing.txt', xd, delimiter=' ')
np.savetxt(filepath + str(filename)+'energy.txt', xd, delimiter=' ')
np.savetxt(filepath + str(filename)+'Q.txt', xQ, delimiter=' ')

def writeList(filename, filelist):
    np.savetxt(filename, filelist, delimiter=" ", fmt = '%s')
 

filelist = glob.glob(filepath+filename+'*ch03*.x_y')
filelist = [elem.split(os.sep)[-1] for elem in filelist]
writeList(filepath + 'filename.txt', filelist)

x=np.zeros([7,2047])
for i in range(len(x)):
    for j in range(2047):
        x[i,j]=j


for i, j in enumerate(channels):
    for jj in range(2047):
        x[i,jj]=Convert_channel_E(x[i,jj], i)


#You can uncomment to create 2D 
"""
dirName = filepath + filename+'image'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ") 
except FileExistsError:
    print("Directory " , dirName ,  " already exists")


for ii in range(real_filenumber):
    figure(ii+3, figsize = [6.4, 1.9])
    imshow(DATA[:,:,ii],cmap='hot', vmin = 0, vmax = DATA.max(), extent = [15,90, 0, 6])
    # xlabel('d-spacing (Angstrom)')
    xlabel('Energy (keV)')
    ylabel('Channel')
    axis('tight')
    if ii<10:
        savefig(dirName+'/'+filename+'_Image_%s.tiff' %ii)   
    else:
        savefig(dirName+'/'+filename+'_Image_%s.tiff' %ii)   
"""
 
