import astropy, photutils, glob, ccdproc,scipy
from astropy.io import fits
from astropy  import units as u
import numpy as np
from astropy.stats import sigma_clipped_stats, SigmaClip, mad_std
from matplotlib import pyplot as plt
from astropy.wcs import WCS
from ccdproc import wcs_project

def binning(file):
    if type(file) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        pass
    else :
        raise Exception("file is not a ccd data")
    file = np.array(file)
    file = np.array(file,dtype = 'float64')
    new_file = []
    for j in range(len(file)):
        row_file = []
        for i in range(1,len(file[j]),2):
            row_file.append(np.mean(file[j][i-1:i+2]))
        row_file = np.array(row_file)
        new_file.append(row_file)
    new_file = np.array(new_file)
    file = new_file.transpose()
    new_file = []
    for j in range(len(file)):
        row_file = []
        for i in range(1,len(file[j]),2):
            row_file.append(np.mean(file[j][i-1:i+2]))
        row_file = np.array(row_file)
        new_file.append(row_file)
    new_file = np.array(new_file)
    new_file = new_file.transpose()
    new_file = ccdproc.CCDData(new_file,unit='adu')
    return new_file

def masterbias(biaslist):
    if type(biaslist) == type([2]):
        pass
    else :
        raise Exception("bias list is not a list")
    for i in biaslist:
            if type(i) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
                pass
            else :
                string = "this list component is not cccd data  "+str(i)

    masterbias = ccdproc.combine(biaslist,method='average',sigma_clip=True,sigma_clip_low_thresh=5,sigma_clip_high_thresh=5,sigma_clip_func=np.ma.median,sigma_clip_dev_func=mad_std)
    return masterbias

def masterdark(darklist,masterbias=None):
    if type(darklist) == type([2]):
        pass
    else :
        raise Exception("dark list is not a list")
    for i in darklist:
            if type(i) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
                pass
            else :
                string = "this list component is not cccd data  "+str(i)

    if  masterbias != None and type(masterbias) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        pass
    elif masterbias != None and type(masterbias) != type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        raise Exception("masterbias is not CCD data")
 
    dlist=[]
    for i in darklist:
      data = i
      if masterbias == None:
           pass
      elif masterbias != None:
          data=ccdproc.subtract_bias(data,masterbias)
      header=data.header
      dlist.append(data)

    masterdark = ccdproc.combine(dlist,method='average',sigma_clip=True,sigma_clip_low_thresh=5,sigma_clip_high_thresh=5,sigma_clip_func=np.ma.median,sigma_clip_dev_func=mad_std)
    header = darklist[0].header
    masterdark.header = header
    return masterdark

def masterflat(flatlist,masterbias=None,masterdark=None):
    if type(flatlist) ==  type([2]):
        pass
    else :
        raise Exception("flat list is not a list")
    for i in flatlist:
            if type(i) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
                pass
            else :
                string = "this list component is not cccd data  "+str(i)

    if  masterbias != None and type(masterbias) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        pass
    elif masterbias != None and type(masterbias) != type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        raise Exception("masterbias is not CCD data")
    if  masterdark != None and type(masterdark) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        pass
    elif masterdark != None and type(masterdark) != type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        raise Exception("masterdark is not CCD data")

    flist = []
    for j in flatlist:
        flat = j
        if masterbias == None and masterdark == None:
          pass 
        elif masterbias != None and masterdark == None:
          flat=ccdproc.subtract_bias(flat,masterbias)
        elif masterbias == None and masterdark != None:
          flat = ccdproc.subtract_dark(flat, masterdark,
                                            exposure_time='EXPTIME',
                                            exposure_unit=u.second,
                                            scale=True)
        elif masterbias != None and masterdark != None:
          flat=ccdproc.subtract_bias(flat,masterbias)
          flat = ccdproc.subtract_dark(flat, masterdark,
                                            exposure_time='EXPTIME',
                                            exposure_unit=u.second,
                                            scale=True)
        flist.append(flat)

        
    def inv_median(a):
      return 1/np.median(a)

    masterflat =ccdproc.combine(flist,method='average',scale=inv_median,
                                sigma_clip=True,sigma_clip_low_thresh=5,sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median,sigma_clip_dev_func=mad_std)
    header = flatlist[0].header
    masterflat.header = header
    return masterflat

def science_frame(sc_list,masterbias=None,masterdark=None,masterflat=None):
    if type(sc_list) ==  type([2]):
        pass
    else :
        raise Exception("sc_list is not a list")
    for i in sc_list:
            if type(i) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
                pass
            else :
                string = "this list component is not cccd data  "+str(i)
    if  masterbias != None and type(masterbias) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        pass
    elif masterbias != None and type(masterbias) != type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        raise Exception("masterbias is not CCD data")
    if  masterdark != None and type(masterdark) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        pass
    elif masterdark != None and type(masterdark) != type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        raise Exception("masterdark is not CCD data")
    if  masterflat != None and type(masterflat) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        pass
    elif masterflat != None and type(masterflat) != type(ccdproc.CCDData(([0,0]),unit = 'adu')):
        raise Exception("masterflat is not CCD data")
    
    #data cleaning begins here with classis logic no more no less
    cleaned_list = []
    for sc in sc_list:
        data = sc
        if masterbias == None and masterdark == None and masterflat == None :
            raise Exception("No calibration file is present how you gonna clean the data")
        elif masterbias == None and masterdark == None and masterflat != None :
            data = ccdproc.flat_correct(data, masterflat)
        elif masterbias == None and masterdark != None and masterflat == None:
            data = ccdproc.subtract_dark(data, masterdark,
                                            exposure_time='EXPTIME',
                                            exposure_unit=u.second,
                                            scale=True)
        elif masterbias == None and masterdark != None and masterflat != None:
            data = ccdproc.subtract_dark(data, masterdark,
                                            exposure_time='EXPTIME',
                                            exposure_unit=u.second,
                                            scale=True)
            data = ccdproc.flat_correct(data, masterflat)
        elif masterbias != None and masterdark == None and masterflat == None :
            data=ccdproc.subtract_bias(data,masterbias)
        elif masterbias != None and masterdark == None and masterflat != None :
            data=ccdproc.subtract_bias(data,masterbias)
            data = ccdproc.flat_correct(data, masterflat)
        elif masterbias != None and masterdark != None and masterflat == None :
            data=ccdproc.subtract_bias(data,masterbias)
            data = ccdproc.subtract_dark(data, masterdark,
                                            exposure_time='EXPTIME',
                                            exposure_unit=u.second,
                                            scale=True)
        elif masterbias != None and masterdark != None and masterflat != None :
            data = ccdproc.subtract_bias(data,masterbias)
            data = ccdproc.subtract_dark(data, masterdark,
                                            exposure_time='EXPTIME',
                                            exposure_unit=u.second,
                                            scale=True)
            data = ccdproc.flat_correct(data, masterflat)
        cleaned_list.append(data)
    return cleaned_list
def stacking_no_wcs(cr_list):
    if type(cr_list) == type([2]):
        pass
    else :
        raise Exception("given list of cleaned files is not a list")
    for i in cr_list:
            if type(i) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
                pass
            else :
                string = "this list component is not cccd data  "+str(i)
                raise Exception(string)
    stacked = ccdproc.combine(cr_list,method='average',sigma_clip=True,sigma_clip_low_thresh=5,sigma_clip_high_thresh=5,sigma_clip_func=np.ma.median,sigma_clip_dev_func=mad_std)
    return stacked

def stacking_wcs(cr_list):
    if type(cr_list) == type([2]):
        pass
    else :
        raise Exception("given list of cleaned files is not a list")
    for i in cr_list:
            if type(i) == type(ccdproc.CCDData(([0,0]),unit = 'adu')):
                pass
            else :
                string = "this list component is not cccd data  "+str(i)
    aligned = []
    for image in cr_list:
        img_aligned = wcs_project(image,cr_list[0].wcs)
        aligned.append(img_aligned)
    image = ccdproc.combine(aligned,method='average',sigma_clip=True,sigma_clip_low_thresh=5,sigma_clip_high_thresh=5,sigma_clip_func=np.ma.median,sigma_clip_dev_func=mad_std)
    return image
