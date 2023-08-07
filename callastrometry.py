#!/usr/local/bin/python

'''
callastrometry - contains functions to handle WCS header information
Requires the following modules: os, docopt, astropy
Contains the following functions: scrubwcsheader, callastrometry

astrometry index files: /opt/astrometry/0.5/data/data.astrometry.net/4200/

Adapted from Natalie Price-Jones.

'''

########################## IMPORT PACKAGES ###########################
import os
import numpy as np
from astropy import wcs
from astropy.io import fits
from scipy.ndimage.interpolation import rotate

########################## DATA LISTS ###########################

# List of file extensions produced by astrometry
keep_fileextensions = ["-indx.png",".new"]
dontkeep_fileextensions = [".axy",".corr","-indx.xyls",".match","-objs.png","-ngc.png",".rdls",".solved",".wcs",".fts"]
# List of header keys corresponding to a WCS solution
keys = ["CTYPE1","CRVAL1","CRPIX1","CDELT1","CROTA1",
        "CTYPE2","CRVAL2","CRPIX2","CDELT2","CROTA2",
        "CD1_1","CD1_2","CD2_1","CD2_2","XPIXSZ","YPIXSZ",
        "EQUINOX","EPOCH","PA","RA","DEC"] 
# keeping OBJRA, OBJDEC; Astrometry.net CRVAL1, CRVAL2 do not correspond to center of image as well as these do
# (used for calling APASS)

# additional DIT header keys that correspond to a Simple Imaging Polynomial (SIP) WCS solution
# (the following are polynomial coefficients)
extrakeys = ["TR1_0","TR1_1","TR1_2","TR1_2","TR1_3","TR1_4","TR1_5","TR1_6","TR1_7",
             "TR1_8","TR1_9","TR1_10","TR1_11","TR1_12","TR1_13","TR1_14",
             "TR2_0","TR2_1","TR2_2","TR2_2","TR2_3","TR2_4","TR2_5","TR2_6","TR2_7",
             "TR2_8","TR2_9","TR2_10","TR2_11","TR2_12","TR2_13","TR2_14",]

def scrubwcsheader(fdir,fname,wcsfolderdir,ind=0):
    '''
    Removes all WCS information from header
    fname:      name of file from which to remove WCS solution
    ind:        index of image in fits file
    Returns nothing explicitly, implicitly saves update file
    '''
    
    # Open file and read header
    light = fits.open(fdir+fname)
    header = light[ind].header
    # obtain RA, DEC header values for future guess
    try:
        CDELT2 = header["CDELT2"]
    except (ValueError,KeyError):
        CDELT2 = -1
    data = light[ind].data
    light.close()
    data_Nup = np.copy(data)
    if CDELT2>=0.0:
        parity = "pos"
        pass
    else:
        # flip image so North is up
        data_Nup = np.flipud(data_Nup)
        parity = "pos"
    # Remove WCS keys
    for key in keys:
        try:
            header.remove(key,ignore_missing=True,remove_all=True)
        except (ValueError,DeprecationWarning):
            continue
    for key in extrakeys:
        try:
            header.remove(key,ignore_missing=True,remove_all=True)
        except (ValueError,DeprecationWarning):
            continue
    # Save updated file
    fits.writeto(wcsfolderdir+fname,data_Nup,header,clobber=True)
    return parity

def callastrometry(fname,parity,fconfig,verifyconfig,generate=True,filekeep=True,new2fts=True):
    '''
    Find WCS solution for an image: the transformation between pixels and R.A./dec
    fname:      name of image file on which to perform astrometry
    generate:   Boolean that specifies whether to force astrometry solution or not
                (kwarg, default = False)
    filekeep:   Boolean that specifies whether to keep auxiliary files from
                astrometry solution (kwarg, default = False)
    Returns nothing explicitly, implicitly saves update file
   '''
    # Read in fits header
    header = fits.getheader(fname)
    # If all WCS entries already in header and generation not forced, done
    if set(keys) < set(header.keys()) and generate == False:
        print 'Astrometry already done'
    # Otherwise, do astrometry
    else:
        # Run astrometry.net
        trimname = fname.split(".fts")[0]
        wcsname = fname.replace(".fts",".wcs")
        
        command = "/opt/astrometry/0.5/bin/solve-field --backend-config {0} --parity {1} --tweak-order 3 --scale-units degwidth --scale-low 0.5 --scale-high 1 --no-fits2fits --use-sextractor --cpulimit 120 --overwrite {2}".format(fconfig,parity,fname)
        verify_command = "/opt/astrometry/0.5/bin/solve-field --backend-config {0} --parity {1} --tweak-order 3 --scale-units degwidth --scale-low 0.5 --scale-high 1 --no-fits2fits --use-sextractor --cpulimit 120 --overwrite --verify {2} {3} --continue".format(verifyconfig,parity,wcsname,fname)
        
        os.system(command)
        os.system(verify_command)
        
        if filekeep==False:
            for ext in dontkeep_fileextensions:
                os.system("rm -f "+trimname+ext)
        if new2fts==True:
            newfname = trimname+"-new.fts"
            oldfname = trimname+".new"
            os.system("mv "+oldfname+" "+newfname)
        else:
            newfname  = trimname+".new"
        '''
        # Read in WCS information
        data,lightwcs = fits.getdata(newfname,header=True)
        w = wcs.WCS(lightwcs)
        # Get x,y pixel scales
        scales = wcs.utils.proj_plane_pixel_scales(w)*3600
        # Update header
        lightwcs["PSCALY"]=scales[0]
        lightwcs["PSCALX"]=scales[1]
        # Save updated file
        fits.writeto(newfname,data,lightwcs,clobber=True)
        '''
