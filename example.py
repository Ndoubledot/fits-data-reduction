import astropy, photutils, glob, ccdproc,scipy
from astropy.io import fits
from astropy  import units as u
import numpy as np
from astropy.stats import sigma_clipped_stats, SigmaClip, mad_std
from matplotlib import pyplot as plt
from astropy.wcs import WCS
from ccdproc import wcs_project
import make
import os
import time
start_time =time.time()
papa_cwd = os.getcwd()
def show(string,image):
	mean, median, std = sigma_clipped_stats(image)
	plt.figure(figsize=(7,7))
	plt.title(string)
	plt.imshow(image,vmin = median - 5*std, vmax = median + 15*std,origin='lower',aspect='equal',cmap='gray')
	plt.colorbar()
	plt.show(block=False)
#------------following example is written using data provided by devanand ullas taken using DOT-------------
#---------------------------thank you devanand so much for the data------------------------------

#we must create list consisting data to put in the functions given in make
#glob is used to create list of file names	
science_frames = sorted(glob.glob('IC310*.fits'))

science_list =[]

for file in science_frames:
        data = ccdproc.CCDData.read(file,unit='adu')
        science_list.append(data)
#now we have science frames as ccd data in the list

#similar thing we do for bias frames, flat frames also with dark frames if it is given  
biasfiles = sorted(glob.glob('bias*.fits'))
biaslist = []
for file in biasfiles:
        data = ccdproc.CCDData.read(file,unit='adu')
        biaslist.append(data)

j_flat = sorted(glob.glob('flats*.fits'))

flat_j = []

for file in j_flat:
        data = ccdproc.CCDData.read(file,unit='adu')
        flat_j.append(data)

#masterbias function takes list consisting list and returns masterbias as ccd data        
masterbias = make.masterbias(biaslist)
#Similary we can make masterdark but we havent as it is not required here
#to make masterdark we would simply put masterdark = make.masterdark(darklist,masterbias)

#masterflat function takes the list, masterbias and masterdak to create calibrated masterflat
#however sometimes calibration files are not given or not needed
#in those times we can skip them and mention None wherever they would be required
#here observations are in J band although script can work for any band
masterflatj = make.masterflat(flat_j,masterbias,None)
#we have put None in spot of masterdark as we do not need it there

science_j = make.science_frame(science_list,masterbias,None,masterflatj)
#science_j is list of ccd data of all cleaned frames, all data is cleaned here
#But we stack all those frames to improve signal to noise ratio
master_j = make.stacking_no_wcs(science_j)
#master_j is single frame made by stacking all frames
#if wcs is mentioned then we can use function stacking_wcs(cr_list)
#if wcs is not mentioned in raw science frame which is the case right now 
#we have used stacking_no_wcs here because raw science frame does not have wcs
#although while using this command we need to make sure that alignment is not need to
#if that be the case then frames should be aligned before stacking
master_j = ccdproc.cosmicray_lacosmic(master_j,readnoise=7.5, sigclip=5,satlevel=65535,niter=20,cleantype='meanmask',gain_apply=True)
#this step is done to remove cosmic ray from the masterframe

show('cleaned science frame',np.array(master_j))
show('raw science frame',np.array(science_list[0]))
#show is funtion written to quickly see the frames
