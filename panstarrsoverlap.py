#!/usr/local/bin/python

#########################################################################################################
#                                                                                                       #
#                                 IMPORT CONFIGURATION FILE                                             #
#                                                                                                       #
#########################################################################################################

import yaml

with open("config.yml","r") as fconfig:
    config_data = yaml.safe_load(fconfig)
with open("panstarrsoverlapconfig.yml","r") as fconfig:
    overlapconfig_data = yaml.safe_load(fconfig)

parentdir             = config_data["parentdir"]
plotdir               = parentdir+"plots/"
passbands             = config_data["passbands"].split(" ")                # list of passbands to be used
calplotdir            = plotdir+"calibration/"                             # directory where calibration plots are stored
calibrateddir         = config_data["calibrateddir"]                       # directory where calibrated images will be stored
finalplotdir          = plotdir+"calibrated/"                              # directory where calibrated plots are stored
sexplotdir_cal        = calplotdir+"sextractor/"                           # subdirectory where SExtractor calibration plots are stored
panstarrsplotdir_cal  = calplotdir+"panstarrs/"                            # subdirectory where Pan-STARRS calibration plots are stored
panstarrsphot         = config_data["panstarrsphot"]                       # Pan-STARRS photometry type: "Ap" (aperture), "PSF", or "Kron"
phottype              = overlapconfig_data["phottype"]                     # SExtractor photometry type; one of ISO, ISOCOR, or AUTO
detect_thresh         = float(overlapconfig_data["detect_thresh"])         # single pixel SNR threshold above RMS background for Source Extractor
snr_thresh            = float(overlapconfig_data["detect_thresh"])         # SNR threshold for photometry measurements
posmatch              = float(overlapconfig_data["posmatch"])              # maximum pointing offset for catalogue matching                [deg]
magmin_r              = float(overlapconfig_data["magmin_r"])              # minimum r Pan-STARRS magnitude for zero-point fit             [mag]
magmax_r              = float(overlapconfig_data["magmax_r"])              # maximum rPan-STARRS magnitude for zero-point fit              [mag]
magmin_i              = float(overlapconfig_data["magmin_i"])              # minimum i Pan-STARRS magnitude for zero-point fit             [mag]
magmax_i              = float(overlapconfig_data["magmax_i"])              # maximum i Pan-STARRS magnitude for zero-point fit             [mag]
magmin_z              = float(overlapconfig_data["magmin_z"])              # minimum z Pan-STARRS magnitude for zero-point fit             [mag]
magmax_z              = float(overlapconfig_data["magmax_z"])              # maximum z Pan-STARRS magnitude for zero-point fit             [mag]
magmin_r_short        = float(overlapconfig_data["magmin_r_short"])        # minimum r Pan-STARRS magnitude for zero-point fit             [mag]
magmax_r_short        = float(overlapconfig_data["magmax_r_short"])        # maximum rPan-STARRS magnitude for zero-point fit              [mag]
flagmax               = float(overlapconfig_data["flagmax"])               # maximum Source Extractor internal flag
clip_sigma            = float(overlapconfig_data["clip_sigma"])            # # of standard deviations for upper and lower clipping limits
zpmin_i               = float(overlapconfig_data["zpmin_i"])               # minimum limit on magnitude zeropoint in the i-band images
zpmax_i               = float(overlapconfig_data["zpmax_i"])               # minimum limit on magnitude zeropoint in the i-band images
zpmin_r               = float(overlapconfig_data["zpmin_r"])               # minimum limit on magnitude zeropoint in the r-band images
zpmax_r               = float(overlapconfig_data["zpmax_r"])               # minimum limit on magnitude zeropoint in the r-band images
zpmin_z               = float(overlapconfig_data["zpmin_z"])               # minimum limit on magnitude zeropoint in the z-band images
zpmax_z               = float(overlapconfig_data["zpmax_z"])               # minimum limit on magnitude zeropoint in the z-band images
zpmin_r_short         = float(overlapconfig_data["zpmin_r_short"])         # minimum limit on magnitude zeropoint in the r-band images
zpmax_r_short         = float(overlapconfig_data["zpmax_r_short"])         # minimum limit on magnitude zeropoint in the r-band images
sexcatdir             = config_data["sexcatdir"]                           # directory where SExtractor catalogue files are stored
sexdir                = config_data["sexdir"]                   # directory where Source Extractor files will be stored

#########################################################################################################
#                                                                                                       #
#                                         IMPORT PACKAGES                                               #
#                                                                                                       #
#########################################################################################################


import os
import numpy as np
from scipy import spatial
from astropy.io import fits
from astropy.io import ascii
from math import log10, floor
from astropy import units as u
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit
from reproject import reproject_interp
from astropy.coordinates import SkyCoord, match_coordinates_sky

import warnings
from scipy.optimize import OptimizeWarning

import pylab
import warnings
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, OptimizeWarning
import matplotlib.patches as mpatches
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=15)

#########################################################################################################
#                                                                                                       #
#                                       IMPORT FUNCTIONS                                                #
#                                                                                                       #
#########################################################################################################

import imp # to force reloading

# filefuncts.py
import filefuncts
imp.reload(filefuncts)
from filefuncts import findfolderdir

# collectzipdata.py
import collectzipdata
imp.reload(collectzipdata)
from collectzipdata import appenditemtodict

#########################################################################################################
#                                                                                                       #
#                                       SOURCE EXTRACTOR                                                #
#                                      CONFIGURATION FILES                                              #
#                                                                                                       #
#########################################################################################################


# default.sex configuration file
sextractor_config = '''# EDITED Default configuration file for SExtractor 2.5.0
# EB 2007-08-27
#
 
#-------------------------------- Catalog ------------------------------------
 
CATALOG_NAME     sex.cat        # name of the output catalog
CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
 
#------------------------------- Extraction ----------------------------------
 
DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
DETECT_MINAREA   16             # minimum number of pixels above threshold
THRESH_TYPE      RELATIVE       # detection is measured relative to the RMS background
DETECT_THRESH    {detect_thresh}  # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
 
FILTER           Y              # apply filter for detection (Y or N)?
FILTER_NAME {filter_name}       # filter
 
DEBLEND_NTHRESH  32             # Number of deblending sub-thresholds
DEBLEND_MINCONT  0.005          # Minimum contrast parameter for deblending
 
CLEAN            Y              # Clean spurious detections? (Y or N)?
CLEAN_PARAM      1.0            # Cleaning efficiency
 
MASK_TYPE        CORRECT        # type of detection MASKing: can be one of
                                # NONE, BLANK or CORRECT
 
#------------------------------ Photometry -----------------------------------
 
SATUR_LEVEL      50000.0        # level (in ADUs) at which saturation arises
 
MAG_ZEROPOINT    0.0            # magnitude zero-point
MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
GAIN             1              # detector gain in e-/ADU
PIXEL_SCALE      1.08           # size of pixel in arcsec (0=use FITS WCS info)
 
#------------------------- Star/Galaxy Separation ----------------------------
 
SEEING_FWHM      4.48           # stellar FWHM in arcsec
STARNNW_NAME     default.nnw    # Neural-Network-Weight table filename
 
#------------------------------ Background -----------------------------------
 
BACK_SIZE        100            # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>
 
BACKPHOTO_TYPE   LOCAL          # can be GLOBAL or LOCAL
BACKPHOTO_THICK 30              # thickness of the background LOCAL annulus (*)
BACK_TYPE       AUTO
BACK_VALUE      0
 
#------------------------------ Check Image ----------------------------------
 
CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                # or APERTURES
CHECKIMAGE_NAME  check.fits     # Filename for the check-image
 
#--------------------- Memory (change with caution!) -------------------------
 
MEMORY_OBJSTACK  30000          # number of objects in stack
MEMORY_PIXSTACK  5000000        # number of pixels in stack
MEMORY_BUFSIZE   1024           # number of lines in buffer
 
#----------------------------- Miscellaneous ---------------------------------
 
VERBOSE_TYPE     NORMAL         # can be QUIET, NORMAL or FULL
WRITE_XML        N              # Write XML file (Y/N)?
XML_NAME         sex.xml        # Filename for XML output
'''

# default.param configuration file
sextractor_params = '''NUMBER
ALPHA_J2000
X_WORLD
ERRX2_WORLD
DELTA_J2000
Y_WORLD
ERRY2_WORLD
MAG_ISO
MAGERR_ISO
FLUX_ISO
FLUXERR_ISO
MAG_ISOCOR
MAGERR_ISOCOR
FLUX_ISOCOR
FLUXERR_ISOCOR
MAG_AUTO
MAGERR_AUTO
FLUX_AUTO
FLUXERR_AUTO
FWHM_IMAGE
FWHM_WORLD
ELLIPTICITY
ELONGATION
BACKGROUND
FLAGS
CLASS_STAR
'''

# default.conv configuation file
default_conv = '''CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
'''

# default.nnw configuration file
default_nnw = '''NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:   9 for profile parameters + 1 for seeing.
# outputs:  ``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.
 3 10 10  1
-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00
-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00
-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00
 0.00000e+00 
 1.00000e+00 
'''

#########################################################################################################
#                                                                                                       #
#                                       DEFINE FUNCTIONS                                                #
#                                                                                                       #
#########################################################################################################

def writesexfiles(dir,VERBOSE,detect_thresh=detect_thresh):
    '''
    Creates configuration (config), parameters (param), filter
    (conv), and Neural Network (nnw) files for Source Extractor.
    '''
    
    sextractor_config_name = "default.sex"
    params_name            = "default.param"
    nnw_name               = "default.nnw"
    conv_name              = "default.conv"
    catalog_name           = "default.cat"

    if VERBOSE:
        verbose_type = "NORMAL"
    else:
        verbose_type = "QUIET"
    
    fp = open(dir+sextractor_config_name, "w+")
    fp.write(sextractor_config.format(detect_thresh=detect_thresh,filter_name=conv_name,
        parameters_name=params_name,starnnw_name=nnw_name,verbose_type=verbose_type))
    fp.close()
    
    fp = open(dir+params_name, "w+")
    fp.write(sextractor_params)
    fp.close()
    if VERBOSE:
        print "wrote file "+dir+params_name
    
    fp = open(dir+conv_name, "w+")
    fp.write(default_conv)
    fp.close()
    if VERBOSE:
        print "wrote file "+dir+conv_name
    
    fp = open(dir+nnw_name, "w+")
    fp.write(default_nnw)
    fp.close()
    if VERBOSE:
        print "wrote file "+dir+nnw_name
    
    return sextractor_config_name,params_name,nnw_name,conv_name,catalog_name

def sexcall(f,folder,fpath,objdir,catdir,bgdir,sexbool,CHECKZP=False):
    '''
    Run source extractor on f, producing a catalogue as well as object and background maps
    f       : file name to run source extractor on
    fpath   : directory path to f
    objdir  : directory to save object map in
    catdir  : directory to save catalogue in
    bgdir   : directory to save background map in
    Returns location of output catalog
    '''

    # Split to make file name for catalogue, object map and background map filenames
    fname = f.split(".fts")[0]
    if sexbool==True:
        if CHECKZP==False:
            # Construct source extractor calls
            objsexcall = "sex -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME "+objdir+folder+"/"+fname+"_objects.fts "+fpath+f
            baksexcall = "sex -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME "+bgdir+folder+"/"+fname+"_background.fts "+fpath+f
            os.system(objsexcall)
            os.system(baksexcall)
        elif CHECKZP==True:
            # extract zero-point for internal calibration
            header        = fits.getheader(fpath+f)
            MAG_ZEROPOINT = header["ZP"]
            # Construct source extractor calls
            objsexcall = "sex -MAG_ZEROPOINT "+str(MAG_ZEROPOINT)+" -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME "+objdir+folder+"/"+fname+"_objects.fts "+fpath+f
            baksexcall = "sex -MAG_ZEROPOINT "+str(MAG_ZEROPOINT)+" -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME "+bgdir+folder+"/"+fname+"_background.fts "+fpath+f
            os.system(objsexcall)
            os.system(baksexcall)
    
    return catdir+folder+"/"+fname+".cat"

def fmagdiff(mag_sex,magerr_sex,mag_ps,magerr_ps):
    '''
    Calculates (Pan-STARRS mag - SExtractor mag) and associated uncertainties.
    '''
    
    magdiff        = np.array(mag_ps) - np.array(mag_sex)
    magdiff_error  = np.sqrt((magerr_ps)**2. + (magerr_sex)**2.)
    
    return magdiff,magdiff_error

def magerrtosnr(magerr):
    '''
    Converts magnitude uncertainty to signal-to-noise ratio.
    '''
    magsnr = 1./magerr
    return magsnr

def magcalibration(file,folder,folderdir,panstarrsdata,sexcat,VERBOSE,GENERATE,CHECKZP):
    '''
    Builds KDTree of Pan-STARRS catalog to be compared with Source Extractor catalogue.
    '''
    
    # load SExtractor catalog data
    sexdata      = ascii.read(sexcat)
    x_world      = np.array(sexdata["X_WORLD"])
    errx2_world  = np.array(sexdata["ERRX2_WORLD"])
    y_world      = np.array(sexdata["Y_WORLD"])
    erry2_world  = np.array(sexdata["ERRY2_WORLD"])
    
    # extract SExtractor R.A. and Dec.
    ra_sex        = x_world
    dec_sex       = y_world
    # extract Pan-STARRS R.A. and Dec.
    ra_ps      =  panstarrsdata["raMean"] 
    dec_ps     =  panstarrsdata["decMean"]

    # search Pan-STARRS to find nearest star to each Source Extractor star
    sex_coord                  = SkyCoord(ra=ra_sex*u.degree,  dec=dec_sex*u.degree)
    PS_coord                   = SkyCoord(ra=ra_ps*u.degree,dec=dec_ps*u.degree)
    idx_closest,d2d_closest,_  = match_coordinates_sky(sex_coord,PS_coord,nthneighbor=1)
    d2d_arcsec_closest         = d2d_closest.arcsecond

    if VERBOSE:
        print len(PS_coord),  "stars extracted from Pan-STARRS catalogue"
        print len(sex_coord), "stars found using Source Extractor"
        print len(idx_closest),   "matching stars"
    
    if GENERATE:
        if CHECKZP==False:
            panstarrsplotfolderdir  = panstarrsplotdir_cal+folder+"/"
        elif CHECKZP==True:
            panstarrsplotfolderdir  = calibrateddir+folder+"/"
        # histogram of pointing offsets between Pan-STARRS and SExtractor catalogs
        pointingoffsethist(file,d2d_arcsec_closest,panstarrsplotfolderdir)

    # clean up Pan-STARRS overlap
    ZEROSLOPE,zeropoint,zeropoint_error,GOODZP,CHECKOBSTIME,mag_sex_matches_clean,magerr_sex_matches_clean,mag_ps_matches_clean,magerr_ps_matches_clean,mag_sex_matches,magerr_sex_matches,mag_ps_matches,magerr_ps_matches = fphotmatches(file,folder,folderdir,panstarrsdata,idx_closest,d2d_arcsec_closest,sexcat,panstarrsplotfolderdir,VERBOSE,GENERATE,CHECKZP)
    
    return ZEROSLOPE,zeropoint,zeropoint_error,GOODZP,CHECKOBSTIME,mag_sex_matches_clean,magerr_sex_matches_clean,mag_ps_matches_clean,magerr_ps_matches_clean,mag_sex_matches,magerr_sex_matches,mag_ps_matches,magerr_ps_matches

def fphotmatches(
    file,folder,folderdir,panstarrsdata,idx_closest,d2d_arcsec_closest,sexcat,panstarrsplotfolderdir,VERBOSE,GENERATE,CHECKZP,flagmax=flagmax,panstarrsphot=panstarrsphot):
    '''
    Cleans up Pan-STARRS and Source Extractor photometry data using information stored in the Pan-STARRS overlapconfig.yml configuration file.
    '''

    header            = fits.getheader(folderdir+file)
    passband          = header["FILTER"]
    exptime           = header["EXPTIME"]
    timeobs           = ftimeobs_to24h(header["TIME-OBS"])

    # SExtractor catalog data
    sexdata           = ascii.read(sexcat)
    alpha_j2000       = np.array(sexdata["ALPHA_J2000"])   # Right Ascension (RA) of barycenter               [hh:mm:ss]
    x_world           = np.array(sexdata["X_WORLD"])       # barycenter position along world x-axis           [deg]
    errx2_world       = np.array(sexdata["ERRX2_WORLD"])   # variance of position along world x-axis          [deg**2]
    xerr_world        = np.sqrt(errx2_world)               # RMS uncertainty of position along world x-axis   [deg]
    delta_j2000       = np.array(sexdata["DELTA_J2000"])   # Declination (Dec) of barycenter                  [dd:mm:ss]
    y_world           = np.array(sexdata["Y_WORLD"])       # barycenter position along world y-axis           [deg]
    erry2_world       = np.array(sexdata["ERRY2_WORLD"])   # variance of position along world y-axis          [deg**2]
    yerr_world        = np.sqrt(erry2_world)               # RMS uncertainty of position along world y-axis   [deg]
    mag_iso           = np.array(sexdata["MAG_ISO"])       # isophotal manitude                               [mag]
    magerr_iso        = np.array(sexdata["MAGERR_ISO"])    # RMS uncertainty for ISO magnitude                [mag]
    flux_iso          = np.array(sexdata["FLUX_ISO"])      # flux density of ISO magnitude                    [ADU]
    fluxerr_iso       = np.array(sexdata["FLUX_ISO"])      # RMS uncertainty of ISO flux density              [ADU]
    mag_isocor        = np.array(sexdata["MAG_ISOCOR"])    # corrected isophotal magnitude                    [mag]
    magerr_isocor     = np.array(sexdata["MAGERR_ISOCOR"]) # RMS uncertainty for ISOCOR magnitude             [mag]
    flux_isocor       = np.array(sexdata["FLUX_ISOCOR"])   # flux density of ISOCOR magnitude                 [ADU]
    fluxerr_isocor    = np.array(sexdata["FLUX_ISOCOR"])   # RMS uncertainty of ISOCOR flux density           [ADU]
    mag_auto          = np.array(sexdata["MAG_AUTO"])      # kron-like elliptical aperture magnitude          [mag]
    magerr_auto       = np.array(sexdata["MAGERR_AUTO"])   # RMS uncertainty for AUTO magnitude               [mag]
    flux_auto         = np.array(sexdata["FLUX_AUTO"])     # flux density of AUTO magnitude                   [ADU]
    fluxerr_auto      = np.array(sexdata["FLUX_AUTO"])     # RMS uncertainty of AUTO flux density             [ADU]
    flags             = np.array(sexdata["FLAGS"])         # extraction flags
    class_star        = np.array(sexdata["CLASS_STAR"])    # S/G classification
    
    # obtain Source Extractor photometry type information
    if phottype=="ISO":
        flux        = flux_iso
        fluxerr     = fluxerr_iso
        mag_sex     = mag_iso
        magerr_sex  = magerr_iso
    elif phottype=="ISOCOR":
        flux        = flux_isocor
        fluxerr     = fluxerr_isocor
        mag_sex     = mag_isocor
        magerr_sex  = magerr_isocor
    elif phottype=="AUTO":
        flux        = flux_auto
        fluxerr     = fluxerr_auto
        mag_sex     = mag_auto
        magerr_sex  = magerr_auto
    
    # Pan-STARRS catalogue data
    PS_RA                = panstarrsdata["raMean"]                       # right ascension from single epoch detections (weighted mean) in equinox J2000             [deg]
    PS_RA_err            = panstarrsdata["raMeanErr"]                    # error in right ascension from single epoch detections (weighted mean) in equinox J2000    [arcsec]
    PS_DEC               = panstarrsdata["decMean"]                      # declination from single epoch detections (weighted mean) in equinox J2000                 [deg]
    PS_DEC_err           = panstarrsdata["decMeanErr"]                   # error in declination from single epoch detections (weighted mean) in equinox J2000        [arcsec]
    PS_r_mag             = panstarrsdata["rMean"+panstarrsphot+"Mag"]    # Sloan r' magnitude                                                                        [mag]
    PS_r_magerr          = panstarrsdata["rMean"+panstarrsphot+"MagErr"] # error in Sloan r' magnitude                                                               [mag]
    PS_rFlags            = panstarrsdata["rFlags"]                       # Information flag bitmask for mean object from r filter detections. Values listed in ObjectFilterFlags.
    PS_i_mag             = panstarrsdata["iMean"+panstarrsphot+"Mag"]    # Sloan i' magnitude                                                                        [mag]
    PS_i_magerr          = panstarrsdata["iMean"+panstarrsphot+"MagErr"] # error in Sloan i' magnitude                                                               [mag]
    PS_iFlags            = panstarrsdata["iFlags"]                       # Information flag bitmask for mean object from i filter detections. Values listed in ObjectFilterFlags.
    PS_z_mag             = panstarrsdata["zMean"+panstarrsphot+"Mag"]    # Sloan z' magnitude                                                                        [mag]
    PS_z_magerr          = panstarrsdata["zMean"+panstarrsphot+"MagErr"] # error in Sloan z' magnitude                                                               [mag]
    PS_zFlags            = panstarrsdata["zFlags"]                       # Information flag bitmask for mean object from z filter detections. Values listed in ObjectFilterFlags.
    PS_objInfoFlag       = panstarrsdata["objInfoFlag"]                  # Information flag bitmask indicating details of the photometry. Values listed in ObjectInfoFlags.
    PS_qualityFlag       = panstarrsdata["qualityFlag"]                  # Subset of objInfoFlag denoting whether this object is real or a likely false positive. Values listed in ObjectQualityFlags.
    
    # Pan-STARRS: closest to Source Extractor
    PS_RA_closest                = PS_RA[idx_closest]
    PS_RA_err_closest            = PS_RA_err[idx_closest]
    PS_DEC_closest               = PS_DEC[idx_closest]
    PS_DEC_err_closest           = PS_DEC_err[idx_closest]
    PS_r_mag_closest             = PS_r_mag[idx_closest]
    PS_r_magerr_closest          = PS_r_magerr[idx_closest]
    PS_rFlags_closest            = PS_rFlags[idx_closest]
    PS_i_mag_closest             = PS_i_mag[idx_closest]
    PS_i_magerr_closest          = PS_i_magerr[idx_closest]
    PS_iFlags_closest            = PS_iFlags[idx_closest]
    PS_z_mag_closest             = PS_z_mag[idx_closest]
    PS_z_magerr_closest          = PS_z_magerr[idx_closest]
    PS_zFlags_closest            = PS_zFlags[idx_closest]
    PS_objInfoFlag_closest       = PS_objInfoFlag[idx_closest]
    PS_qualityFlag_closest       = PS_qualityFlag[idx_closest]
    
    # obtain DIT passband information
    if passband=="i":
        mag_ps_closest        = PS_i_mag_closest
        magerr_ps_closest     = PS_i_magerr_closest
        magsnr_ps_closest     = magerrtosnr(PS_i_magerr_closest)
        passband_flag_closest = PS_iFlags_closest
        magmin                = magmin_i
        magmax                = magmax_i
    elif passband=="r":
        mag_ps_closest        = PS_r_mag_closest
        magerr_ps_closest     = PS_r_magerr_closest
        magsnr_ps_closest     = magerrtosnr(PS_r_magerr_closest)
        passband_flag_closest = PS_rFlags_closest
        if exptime==5.:
            magmin                = magmin_r_short
            magmax                = magmax_r_short
        else:
            magmin                = magmin_r
            magmax                = magmax_r
    elif passband=="z":
        mag_ps_closest        = PS_z_mag_closest
        magerr_ps_closest     = PS_z_magerr_closest
        magsnr_ps_closest     = magerrtosnr(PS_z_magerr_closest)
        passband_flag_closest = PS_zFlags_closest
        magmin                = magmin_z
        magmax                = magmax_z

    # find matches with basic quality cleaning
    idx_matches = (d2d_arcsec_closest<=posmatch) & (flux>0.0) & (fluxerr>0.0) & (flags<=3.0) & (mag_ps_closest!=-999.)

    # Pan-STARRS: matches to Source Extractor
    mag_ps_matches       = mag_ps_closest[idx_matches] 
    magerr_ps_matches    = magerr_ps_closest[idx_matches]
    # Source Extractor: matches to Pan-STARRS
    mag_sex_matches      = mag_sex[idx_matches]
    magerr_sex_matches   = magerr_sex[idx_matches]

    # find matches with proper quality cleaning
    idx_matches_clean = (d2d_arcsec_closest<=posmatch) & (flux>0.0) & (fluxerr>0.0) & (flags<=3.0) & (mag_ps_closest<=magmax) & (mag_ps_closest>=magmin) & (magsnr_ps_closest>=snr_thresh) & (passband_flag_closest<8388608.) & (PS_qualityFlag_closest<128.) 

    # Pan-STARRS: matches to Source Extractor
    mag_ps_matches_clean       = mag_ps_closest[idx_matches_clean] 
    magerr_ps_matches_clean    = magerr_ps_closest[idx_matches_clean]
    # Source Extractor: matches to Pan-STARRS
    mag_sex_matches_clean      = mag_sex[idx_matches_clean]
    magerr_sex_matches_clean   = magerr_sex[idx_matches_clean]
    
    # perform sigma-clipping on magnitude zero points
    #ra_sc,dec_sc,x_world_sc,y_world_sc,mag_ps_sc,magerr_ps_sc,mag_sex_sc,magerr_sex_sc,zp_sc=fzeropoint_sc(ra_clean,dec_clean,x_world_clean,y_world_clean,mag_ps_clean,magerr_ps_clean,mag_sex_clean,magerr_sex_clean,clip_sigma)
    
    # perform linear fitting to calculate magnitude zero point and slope
    ZEROSLOPE,zeropoint,zeropoint_error = fitzeropoint(file,passband,exptime,mag_sex_matches_clean,magerr_sex_matches_clean,mag_ps_matches_clean,magerr_ps_matches_clean,mag_sex_matches,magerr_sex_matches,mag_ps_matches,magerr_ps_matches,panstarrsplotfolderdir,CHECKZP=CHECKZP)
    #ZEROSLOPE_sc                        = fitzeropoint(file,passband,exptime,mag_sex_sc,magerr_sex_sc,mag_ps_sc,magerr_ps_sc,mag_sex_justmatches,magerr_sex_justmatches,mag_ps_justmatches,magerr_ps_justmatches,panstarrsplotfolderdir,SC=True)
    
    GOODZP       = fgoodzp(zeropoint,passband,exptime,CHECKZP)
    CHECKOBSTIME = fcheckobstime(passband,folder,timeobs,exptime)
        
    return ZEROSLOPE,zeropoint,zeropoint_error,GOODZP,CHECKOBSTIME,mag_sex_matches_clean,magerr_sex_matches_clean,mag_ps_matches_clean,magerr_ps_matches_clean,mag_sex_matches,magerr_sex_matches,mag_ps_matches,magerr_ps_matches

def pointingoffsethist(file,dist,plotdir,deltapos=1.0,poslower=0.0,posupper=10.0):
    '''
    Plots the pointing offset between Pan-STARRS and Source Extractor catalog matches.
    '''
    
    newfname = file.replace(".fts","_PS_posresidualhist.pdf")
    
    N = len(np.arange(poslower,posupper,deltapos))
    params = dict(bins=N,range=(poslower,posupper))
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    pylab.hist(dist,**params)
    plt.axvline(posmatch,linestyle="dashed",color="black",label="position match")
    plt.xlabel("PS Residual Pointing Offset (arcseconds)")
    plt.ylabel("N")
    plt.legend(loc="upper right")
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)

def fzeropoint_sc(ra,dec,x_world,y_world,mag_ps,magerr_ps,mag_sex,magerr_sex,clip_sigma):
    '''
    Measures the magnitude zero point of an image.
    '''
    
    zp            = mag_sex - mag_ps
    zp_sigmaclip  = sigma_clip(zp,sigma=clip_sigma,iters=5)
    zp_outliers   = zp_sigmaclip.mask
    zp_sc         = zp_sigmaclip[~zp_outliers].data
    
    # clean up using sigma-clipped mask
    ra_sc             = np.ma.masked_array(ra,mask=zp_outliers)[~zp_outliers].data
    dec_sc            = np.ma.masked_array(dec,mask=zp_outliers)[~zp_outliers].data
    x_world_sc        = np.ma.masked_array(x_world,mask=zp_outliers)[~zp_outliers].data
    y_world_sc        = np.ma.masked_array(y_world,mask=zp_outliers)[~zp_outliers].data
    mag_ps_sc         = np.ma.masked_array(mag_ps,mask=zp_outliers)[~zp_outliers].data
    magerr_ps_sc      = np.ma.masked_array(magerr_ps,mask=zp_outliers)[~zp_outliers].data
    mag_sex_sc        = np.ma.masked_array(mag_sex,mask=zp_outliers)[~zp_outliers].data
    magerr_sex_sc     = np.ma.masked_array(magerr_sex,mask=zp_outliers)[~zp_outliers].data

    return ra_sc,dec_sc,x_world_sc,y_world_sc,mag_ps_sc,magerr_ps_sc,mag_sex_sc,magerr_sex_sc,zp_sc

def fitzeropoint(file,passband,exptime,mag_sex_matches_clean,magerr_sex_matches_clean,mag_ps_matches_clean,magerr_ps_matches_clean,mag_sex_matches,magerr_sex_matches,mag_ps_matches,magerr_ps_matches,plotdir,SC=False,YERRORS=False,PLOTMEDIAN=True,PLOTWEIGHTEDAVG=False,CHECKZP=False):
    '''
    Fits a straight line to (SExtractor - Pan-STARRS) vs. Pan-STARRS magnitudes to measure magnitude zero point and slope in the data.
    '''
    
    def fline(mag_ps_matches_clean,slope,zeropoint):
        y =  slope*mag_ps_matches_clean + zeropoint
        return y
    
    def fline_odr(params,mag_ps_matches_clean):
        slope = params[0]
        zeropoint = params[1]
        y = slope*mag_ps_matches_clean + zeropoint
        return y
    
    if passband=="r":
        if exptime==5.:
            magmin = magmin_r_short
            magmax = magmax_r_short
        else:
            magmin = magmin_r
            magmax = magmax_r
    elif passband=="i":
        magmin = magmin_i
        magmax = magmax_i
    elif passband=="z":
        magmin = magmin_z
        magmax = magmax_z
    
    # raw data (for plotting)
    magdiff_matches,magdiff_error_matches = fmagdiff(mag_sex_matches,magerr_sex_matches,mag_ps_matches,magerr_ps_matches)
    
    # only fit data if there's more than two stars
    if len(mag_sex_matches_clean)>2:
        magdiff,magdiff_error       = fmagdiff(mag_sex_matches_clean,magerr_sex_matches_clean,mag_ps_matches_clean,magerr_ps_matches_clean)
        magdiff_median              = np.median(magdiff)
        magdiff_median_error        = np.std(magdiff)/np.sqrt(len(magdiff))
        magdiff_weights             = 1./(magdiff_error)**2.
        magdiff_weightedavg         = np.average(magdiff,weights=magdiff_weights)
        magdiff_median_rounded,magdiff_median_error_rounded = round_with_uncertainties(magdiff_median,magdiff_median_error)
    
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            try:
                CURVEFIT=True
                # ZP should be zero if checking calibrated frames
                if CHECKZP==False:
                    guess = np.array([0.0,25.5])
                elif CHECKZP==True:
                    guess = np.array([0.0,0.0])

                if YERRORS==False:
                    popt,pcov           = curve_fit(fline,mag_ps_matches_clean,magdiff,p0=guess,maxfev=1000)
                    popt_uncertainties  = np.sqrt(np.diag(pcov))
                elif YERRORS==True:
                    popt,pcov            = curve_fit(fline,mag_ps_matches_clean,magdiff,p0=guess,sigma=magdiff_error,absolute_sigma=True,maxfev=1000)
                    popt_uncertainties   = np.sqrt(np.diag(pcov))
                maglower                 = np.min(mag_ps_matches_clean)
                magupper                 = np.max(mag_ps_matches_clean)
                mag_ps_fit               = np.linspace(maglower,magupper,100)
                deltamag_fit             = fline(mag_ps_fit,*popt)
                
                slope                    = popt[0]
                slope_error              = popt_uncertainties[0]
                zeropoint                = popt[1]
                zeropoint_error          = popt_uncertainties[1]
                
                slope_rounded,slope_error_rounded                   = round_with_uncertainties(slope,slope_error)
                zeropoint_rounded,zeropoint_error_rounded           = round_with_uncertainties(zeropoint,zeropoint_error)
                
                errorsig = 3.0
                ZEROSLOPE = isclose(0.0,slope,abs_tol=errorsig*slope_error)
        
            except (RuntimeError, OptimizeWarning):
                CURVEFIT  = False
                ZEROSLOPE = False
    else:
        magdiff_median       = None
        magdiff_median_error = None
        CURVEFIT             = False
        ZEROSLOPE            = False
    
    # scatterplot
    if SC==True:
        newfname = file.replace(".fts","_PS_magsoln_zeropoints_sc.pdf")
    else:
        newfname = file.replace(".fts","_PS_magsoln_zeropoints.pdf")
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    # raw data 
    plt.errorbar(mag_ps_matches,magdiff_matches,xerr=magerr_ps_matches,yerr=magdiff_error_matches,fmt="o",color="blue",alpha=0.5,label="data (raw)")
    # cleaned data
    if len(mag_sex_matches_clean)>2:
        plt.errorbar(mag_ps_matches_clean,magdiff,xerr=magerr_ps_matches_clean,yerr=magdiff_error,fmt="o",color="red",alpha=0.5,label="data (clean)")
        if PLOTMEDIAN==True:
            plt.axhline(magdiff_median,color="black",linestyle="--",linewidth=2,label="median (clean)")
        plt.text(0.05,0.85,"zeropoint: "+str(magdiff_median_rounded)+" +/- "+str(magdiff_median_error_rounded),transform=ax.transAxes)
        if PLOTWEIGHTEDAVG==True:
            plt.axhline(magdiff_weightedavg,color="black",linestyle="-.",linewidth=2,label="weighted average (clean)")
    if CURVEFIT==True:
        if YERRORS==True:
            plt.plot(mag_ps_fit,deltamag_fit,color="purple",linewidth=3,label="curvefit with y errors")
        elif YERRORS==False:
            plt.plot(mag_ps_fit,deltamag_fit,color="purple",linewidth=3,label="curvefit without y errors")
        plt.text(0.05,0.95,"slope: "+str(slope_rounded)+" +/- "+str(slope_error_rounded),transform=ax.transAxes)
        if ZEROSLOPE==True:
            plt.text(0.05,0.9,"slope is zero",transform=ax.transAxes)
    plt.axvline(magmin,linestyle="dashed",color="black")
    plt.axvline(magmax,linestyle="dashed",color="black")
    plt.xlabel("Pan-STARRS Magnitude")
    plt.ylabel("(Pan-STARRS - SExtractor Magnitude)")
    plt.legend(prop={"size":10},loc="lower left",numpoints=1)
    plt.xlim(8.,22.)
    if CHECKZP==False:
        plt.ylim(20.,34.)
    elif CHECKZP==True:
        plt.ylim(-10.,10.)
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)
    
    if CURVEFIT==True:
        print "slope: "+str(slope_rounded)+" +/- "+str(slope_error_rounded)
        if ZEROSLOPE==True:
            print "slope is zero to within uncertainties"
        else:
            print "slope is not zero to wihin uncertainties"
        print "zeropoint: "+str(magdiff_median_rounded)+" +/- "+str(magdiff_median_error_rounded)
    else:
        print "unsuccessful fit"
    
    return ZEROSLOPE,magdiff_median,magdiff_median_error

def plotzp_all(mag_sex_matches_clean_dict,magerr_sex_matches_clean_dict,mag_ps_matches_clean_dict,magerr_ps_matches_clean_dict,mag_sex_matches_dict,magerr_sex_matches_dict,mag_ps_matches_dict,magerr_ps_matches_dict,plotdir,exp,MULTIPROCESSING_CAL,CHECKZP):
    # r-band
    if "r" in mag_sex_matches_dict.keys():
        mag_sex_r_matches                                     = np.array(mag_sex_matches_dict["r"])
        magerr_sex_r_matches                                  = np.array(magerr_sex_matches_dict["r"])
        mag_ps_r_matches                                      = np.array(mag_ps_matches_dict["r"])
        magerr_ps_r_matches                                   = np.array(magerr_ps_matches_dict["r"])
        magdiff_r_matches,magdiff_err_r_matches               = fmagdiff(mag_sex_r_matches,magerr_sex_r_matches,mag_ps_r_matches,magerr_ps_r_matches)
    if "r" in mag_sex_matches_clean_dict.keys():
        mag_sex_r_matches_clean                               = np.array(mag_sex_matches_clean_dict["r"])
        magerr_sex_r_matches_clean                            = np.array(magerr_sex_matches_clean_dict["r"])
        mag_ps_r_matches_clean                                = np.array(mag_ps_matches_clean_dict["r"])
        magerr_ps_r_matches_clean                             = np.array(magerr_ps_matches_clean_dict["r"])
        magdiff_r_matches_clean,magdiff_err_r_matches_clean   = fmagdiff(mag_sex_r_matches_clean,magerr_sex_r_matches_clean,mag_ps_r_matches_clean,magerr_ps_r_matches_clean)
        magdiff_median_r_matches_clean                        = np.median(magdiff_r_matches_clean)
        magdiff_err_median_r_matches_clean                    = np.std(magdiff_r_matches_clean)/np.sqrt(len(magdiff_r_matches_clean))
    
    # i-band
    if "i" in mag_sex_matches_dict.keys():
        mag_sex_i_matches                                     = np.array(mag_sex_matches_dict["i"])
        magerr_sex_i_matches                                  = np.array(magerr_sex_matches_dict["i"])
        mag_ps_i_matches                                      = np.array(mag_ps_matches_dict["i"])
        magerr_ps_i_matches                                   = np.array(magerr_ps_matches_dict["i"])
        magdiff_i_matches,magdiff_err_i_matches               = fmagdiff(mag_sex_i_matches,magerr_sex_i_matches,mag_ps_i_matches,magerr_ps_i_matches)
    if "i" in mag_sex_matches_clean_dict.keys():
        mag_sex_i_matches_clean                               = np.array(mag_sex_matches_clean_dict["i"])
        magerr_sex_i_matches_clean                            = np.array(magerr_sex_matches_clean_dict["i"])
        mag_ps_i_matches_clean                                = np.array(mag_ps_matches_clean_dict["i"])
        magerr_ps_i_matches_clean                             = np.array(magerr_ps_matches_clean_dict["i"])
        magdiff_i_matches_clean,magdiff_err_i_matches_clean   = fmagdiff(mag_sex_i_matches_clean,magerr_sex_i_matches_clean,mag_ps_i_matches_clean,magerr_ps_i_matches_clean)
        magdiff_median_i_matches_clean                        = np.median(magdiff_i_matches_clean)
        magdiff_err_median_i_matches_clean                    = np.std(magdiff_i_matches_clean)/np.std(len(magdiff_i_matches_clean))

    # z-band
    if "z" in mag_sex_matches_dict.keys():
        mag_sex_z_matches                                     = np.array(mag_sex_matches_dict["z"])
        magerr_sex_z_matches                                  = np.array(magerr_sex_matches_dict["z"])
        mag_ps_z_matches                                      = np.array(mag_ps_matches_dict["z"])
        magerr_ps_z_matches                                   = np.array(magerr_ps_matches_dict["z"])
        magdiff_z_matches,magdiff_err_z_matches               = fmagdiff(mag_sex_z_matches,magerr_sex_z_matches,mag_ps_z_matches,magerr_ps_z_matches)
    if "z" in mag_sex_matches_clean_dict.keys():
        mag_sex_z_matches_clean                               = np.array(mag_sex_matches_clean_dict["z"])
        magerr_sex_z_matches_clean                            = np.array(magerr_sex_matches_clean_dict["z"])
        mag_ps_z_matches_clean                                = np.array(mag_ps_matches_clean_dict["z"])
        magerr_ps_z_matches_clean                             = np.array(magerr_ps_matches_clean_dict["z"])
        magdiff_z_matches_clean,magdiff_err_z_matches_clean   = fmagdiff(mag_sex_z_matches_clean,magerr_sex_z_matches_clean,mag_ps_z_matches_clean,magerr_ps_z_matches_clean)
        magdiff_median_z_matches_clean                        = np.median(magdiff_z_matches_clean)
        magdiff_err_median_z_matches_clean                    = np.std(magdiff_z_matches_clean)/np.sqrt(len(magdiff_z_matches_clean))
    
    if ("r" in mag_sex_matches_dict.keys()) or ("r" in mag_sex_matches_clean_dict.keys()):

        # scatterplot
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # just matches
        if "r" in mag_sex_matches_dict.keys():
            plt.errorbar(mag_ps_r_matches,magdiff_r_matches,xerr=magerr_ps_r_matches,yerr=magdiff_err_r_matches,fmt="o",markersize=8,color="blue",alpha=0.5,label="data (all matches)")
        # clean
        if "r" in mag_sex_matches_clean_dict.keys():
            plt.errorbar(mag_ps_r_matches_clean,magdiff_r_matches_clean,xerr=magerr_ps_r_matches_clean,yerr=magdiff_err_r_matches_clean,fmt="o",markersize=8,color="red",alpha=0.5,label="data (clean)")
            plt.axhline(magdiff_median_r_matches_clean,color="black",linestyle="--",linewidth=2,label="median")
            plt.text(0.05,0.85,"ZP: "+str(magdiff_median_r_matches_clean)+" +/- "+str(round(magdiff_err_median_r_matches_clean,3))+" (median)",transform=ax.transAxes)
        plt.xlabel("Pan-STARRS Magnitude")
        plt.ylabel("(SExtractor Magnitude - Pan-STARRS Magnitude)")
        plt.legend(prop={"size":10},loc="lower left",numpoints=1)
        plt.xlim(8,22)
        if CHECKZP==False:
            plt.ylim(20,34)
        elif CHECKZP==True:
            plt.ylim(-10.,10.)
        if exp=="SHORT":
            plt.axvline(magmin_r_short,linestyle="dashed",color="black")
            plt.axvline(magmax_r_short,linestyle="dashed",color="black")
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_all_r_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_all_r_shortexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="LONG":
            plt.axvline(magmin_r,linestyle="dashed",color="black")
            plt.axvline(magmax_r,linestyle="dashed",color="black")
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_all_r_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_all_r_longexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # 2d histogram
        # clean
        if "r" in mag_sex_matches_clean_dict.keys():
            fig = plt.figure(1,figsize=(11,8.5))
            ax = fig.add_subplot(111)
            plt.hist2d(mag_ps_r_matches_clean,magdiff_r_matches_clean,bins=(40,20),cmin=1,cmap=plt.cm.plasma)
            plt.axhline(magdiff_median_r_matches_clean,color="black",linestyle="--",linewidth=2,label="median")
            plt.text(0.05,0.85,"ZP: "+str(magdiff_median_r_matches_clean)+" +/- "+str(round(magdiff_err_median_r_matches_clean,3))+" (median)",transform=ax.transAxes)
            plt.xlabel("Pan-STARRS Magnitude")
            plt.ylabel("(SExtractor Magnitude - Pan-STARRS Magnitude)")
            plt.legend(prop={"size":10},loc="lower left",numpoints=1)
            plt.colorbar()
            plt.xlim(8,22)
            if CHECKZP==False:
                plt.ylim(20,34)
            elif CHECKZP==True:
                plt.ylim(-10.,10.)
            if exp=="SHORT":
                plt.axvline(magmin_r_short,linestyle="dashed",color="black")
                plt.axvline(magmax_r_short,linestyle="dashed",color="black")
                if CHECKZP==False:
                    plt.savefig(plotdir+"PS_zp_all_r_shortexp_2dhist.pdf",bbox_inches="tight")
                elif CHECKZP==True:
                    plt.savefig(plotdir+"PS_zp_all_r_shortexp_2dhist_calibrated.pdf",bbox_inches="tight")
            elif exp=="LONG":
                plt.axvline(magmin_r,linestyle="dashed",color="black")
                plt.axvline(magmax_r,linestyle="dashed",color="black")
                if CHECKZP==False:
                    plt.savefig(plotdir+"PS_zp_all_r_longexp_2dhist.pdf",bbox_inches="tight")
                elif CHECKZP==True:
                    plt.savefig(plotdir+"PS_zp_all_r_longexp_2dhist_calibrated.pdf",bbox_inches="tight")
            plt.close(fig)

    if ("i" in mag_sex_matches_dict.keys()) or ("i" in mag_sex_matches_clean_dict.keys()):

        # scatterplot
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # just matches
        if "i" in mag_sex_matches_dict.keys():
            plt.errorbar(mag_ps_i_matches,magdiff_i_matches,xerr=magerr_ps_i_matches,yerr=magdiff_err_i_matches,fmt="o",markersize=8,color="blue",alpha=0.5,label="data (all matches)")
        # clean
        if "i" in mag_sex_matches_clean_dict.keys():
            plt.errorbar(mag_ps_i_matches_clean,magdiff_i_matches_clean,xerr=magerr_ps_i_matches_clean,yerr=magdiff_err_i_matches_clean,fmt="o",markersize=8,color="red",alpha=0.5,label="data (clean)")
            plt.axhline(magdiff_median_i_matches_clean,color="black",linestyle="--",linewidth=2,label="median")
            plt.text(0.05,0.85,"ZP: "+str(magdiff_median_i_matches_clean)+" +/- "+str(round(magdiff_err_median_i_matches_clean,3))+" (median)",transform=ax.transAxes)
        plt.xlabel("Pan-STARRS Magnitude")
        plt.ylabel("(SExtractor Magnitude - Pan-STARRS Magnitude)")
        plt.legend(prop={"size":10},loc="lower left",numpoints=1)
        plt.axvline(magmin_i,linestyle="dashed",color="black")
        plt.axvline(magmax_i,linestyle="dashed",color="black")
        plt.xlim(8,22)
        if CHECKZP==False:
            plt.ylim(20,34)
        elif CHECKZP==True:
            plt.ylim(-10.,10.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_all_i_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_all_i_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_all_i_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_all_i_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # 2d histogram
        # clean
        if "i" in mag_sex_matches_clean_dict.keys():
            fig = plt.figure(1,figsize=(11,8.5))
            ax = fig.add_subplot(111)
            plt.hist2d(mag_ps_i_matches_clean,magdiff_i_matches_clean,bins=(40,20),cmin=1,cmap=plt.cm.plasma)
            plt.axhline(magdiff_median_i_matches_clean,color="black",linestyle="--",linewidth=2,label="median")
            plt.text(0.05,0.85,"ZP: "+str(magdiff_median_i_matches_clean)+" +/- "+str(round(magdiff_err_median_i_matches_clean,3))+" (median)",transform=ax.transAxes)
            plt.xlabel("Pan-STARRS Magnitude")
            plt.ylabel("(SExtractor Magnitude - Pan-STARRS Magnitude)")
            plt.legend(prop={"size":10},loc="lower left",numpoints=1)
            plt.axvline(magmin_i,linestyle="dashed",color="black")
            plt.axvline(magmax_i,linestyle="dashed",color="black")
            plt.colorbar()
            plt.xlim(8,22)
            if CHECKZP==False:
                plt.ylim(20,34)
            elif CHECKZP==True:
                plt.ylim(-10.,10.)
            if exp=="LONG":
                if CHECKZP==False:
                    plt.savefig(plotdir+"PS_zp_all_i_longexp_2dhist.pdf",bbox_inches="tight")
                elif CHECKZP==True:
                    plt.savefig(plotdir+"PS_zp_all_i_longexp_2dhist_calibrated.pdf",bbox_inches="tight")
            elif exp=="SHORT":
                if CHECKZP==False:
                    plt.savefig(plotdir+"PS_zp_all_i_shortexp_2dhist.pdf",bbox_inches="tight")
                elif CHECKZP==True:
                    plt.savefig(plotdir+"PS_zp_all_i_shortexp_2dhist_calibrated.pdf",bbox_inches="tight")
            plt.close(fig)
    
    if ("z" in mag_sex_matches_dict.keys()) or ("z" in mag_sex_matches_clean_dict.keys()):

        # scatterplot
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # just matches
        if "z" in mag_sex_matches_dict.keys():
            plt.errorbar(mag_ps_z_matches,magdiff_z_matches,xerr=magerr_ps_z_matches,yerr=magdiff_err_z_matches,fmt="o",markersize=8,color="blue",alpha=0.5,label="data (all matches)")
        # clean
        if "z" in mag_sex_matches_clean_dict.keys():
            plt.errorbar(mag_ps_z_matches_clean,magdiff_z_matches_clean,xerr=magerr_ps_z_matches_clean,yerr=magdiff_err_z_matches_clean,fmt="o",markersize=8,color="red",alpha=0.5,label="data (clean)")
            plt.axhline(magdiff_median_z_matches_clean,color="black",linestyle="--",linewidth=2,label="median")
            plt.text(0.05,0.85,"ZP: "+str(magdiff_median_z_matches_clean)+" +/- "+str(round(magdiff_err_median_z_matches_clean,3))+" (median)",transform=ax.transAxes)
        plt.xlabel("Pan-STARRS Magnitude")
        plt.ylabel("(SExtractor Magnitude - Pan-STARRS Magnitude)")
        plt.legend(prop={"size":10},loc="lower left",numpoints=1)
        plt.axvline(magmin_z,linestyle="dashed",color="black")
        plt.axvline(magmax_z,linestyle="dashed",color="black")
        plt.xlim(8,22)
        if CHECKZP==False:
            plt.ylim(20,34)
        elif CHECKZP==True:
            plt.ylim(-10.,10.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_all_z_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_all_z_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_all_z_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_all_z_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # 2d histogram
        # clean
        if "z" in mag_sex_matches_clean_dict.keys():
            fig = plt.figure(1,figsize=(11,8.5))
            ax = fig.add_subplot(111)
            plt.hist2d(mag_ps_z_matches_clean,magdiff_z_matches_clean,bins=(40,20),cmin=1,cmap=plt.cm.plasma)
            plt.axhline(magdiff_median_z_matches_clean,color="black",linestyle="--",linewidth=2,label="median")
            plt.text(0.05,0.85,"ZP: "+str(magdiff_median_z_matches_clean)+" +/- "+str(round(magdiff_err_median_z_matches_clean,3))+" (median)",transform=ax.transAxes)
            plt.xlabel("Pan-STARRS Magnitude")
            plt.ylabel("(SExtractor Magnitude - Pan-STARRS Magnitude)")
            plt.legend(prop={"size":10},loc="lower left",numpoints=1)
            plt.axvline(magmin_z,linestyle="dashed",color="black")
            plt.axvline(magmax_z,linestyle="dashed",color="black")
            plt.colorbar()
            plt.xlim(8,22)
            if CHECKZP==False:
                plt.ylim(20,34)
            elif CHECKZP==True:
                plt.ylim(-10.,10.)
            if exp=="LONG":
                if CHECKZP==False:
                    plt.savefig(plotdir+"PS_zp_all_z_longexp_2dhist.pdf",bbox_inches="tight")
                elif CHECKZP==True:
                    plt.savefig(plotdir+"PS_zp_all_z_longexp_2dhist_calibrated.pdf",bbox_inches="tight")
            elif exp=="SHORT":
                if CHECKZP==False:
                    plt.savefig(plotdir+"PS_zp_all_z_shortexp_2dhist.pdf",bbox_inches="tight")
                elif CHECKZP==True:
                    plt.savefig(plotdir+"PS_zp_all_z_shortexp_2dhist_calibrated.pdf",bbox_inches="tight")
            plt.close(fig)

def plotzpcorrs(zp_dict,zp_err_dict,zp_all_dict,zp_err_all_dict,airmass_dict,airmass_all_dict,timeobs_dict,timeobs_all_dict,FWHM_dict,FWHM_all_dict,exp,plotdir,CHECKZP):
    '''
    Plots correlation between zero-point and airmass.
    '''
    
    # r-band
    if "r" in zp_dict.keys():
        zp_r          = np.array(zp_dict["r"])
        zp_err_r      = np.array(zp_err_dict["r"])
        airmass_r     = np.array(airmass_dict["r"])
        timeobs_r     = np.array(timeobs_dict["r"])
        FWHM_r        = np.array(FWHM_dict["r"])
    if "r" in zp_all_dict.keys():
        zp_all_r      = np.array(zp_all_dict["r"])
        zp_err_all_r  = np.array(zp_err_all_dict["r"])
        airmass_all_r = np.array(airmass_all_dict["r"])
        timeobs_all_r = np.array(timeobs_all_dict["r"])
        FWHM_all_r    = np.array(FWHM_all_dict["r"])

    # i-band
    if "i" in zp_dict.keys():
        zp_i          = np.array(zp_dict["i"])
        zp_err_i      = np.array(zp_err_dict["i"])
        airmass_i     = np.array(airmass_dict["i"])
        timeobs_i     = np.array(timeobs_dict["i"])
        FWHM_i        = np.array(FWHM_dict["i"])
    if "i" in zp_all_dict.keys():
        zp_all_i      = np.array(zp_all_dict["i"])
        zp_err_all_i  = np.array(zp_err_all_dict["i"])
        airmass_all_i = np.array(airmass_all_dict["i"])
        timeobs_all_i = np.array(timeobs_all_dict["i"])
        FWHM_all_i    = np.array(FWHM_all_dict["i"])

    # z-band
    if "z" in zp_dict.keys():
        zp_z          = np.array(zp_dict["z"])
        zp_err_z      = np.array(zp_err_dict["z"])
        airmass_z     = np.array(airmass_dict["z"])
        timeobs_z     = np.array(timeobs_dict["z"])
        FWHM_z        = np.array(FWHM_dict["z"])
    if "z" in zp_all_dict.keys():
        zp_all_z      = np.array(zp_all_dict["z"])
        zp_err_all_z  = np.array(zp_err_all_dict["z"])
        airmass_all_z = np.array(airmass_all_dict["z"])
        timeobs_all_z = np.array(timeobs_all_dict["z"])
        FWHM_all_z    = np.array(FWHM_all_dict["z"])
    
    if ("r" in zp_dict.keys()) or ("r" in zp_all_dict.keys()):

        # airmass correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "r" in zp_all_dict.keys():
            plt.errorbar(airmass_all_r,zp_all_r,yerr=zp_err_all_r,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "r" in zp_dict.keys():
            plt.errorbar(airmass_r,zp_r,yerr=zp_err_r,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Airmass (multiple of zenithal airmass)")
        plt.ylabel("PS Zero-point (mags)")
        plt.xlim(1.1,1.35)
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_airmass_r_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_airmass_r_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_airmass_r_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_airmass_r_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # timeobs correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "r" in zp_all_dict.keys():
            plt.errorbar(timeobs_all_r,zp_all_r,yerr=zp_err_all_r,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "r" in zp_dict.keys():
            plt.errorbar(timeobs_r,zp_r,yerr=zp_err_r,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Time of Observtion (hh.hh)")
        plt.ylabel("PS Zero-point (mags)")
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_r_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_r_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_r_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_r_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # FWHM correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "r" in zp_all_dict.keys():
            plt.errorbar(FWHM_all_r,zp_all_r,yerr=zp_err_all_r,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "r" in zp_dict.keys():
            plt.errorbar(FWHM_r,zp_r,yerr=zp_err_r,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Mean PSF FWHM (pixels)")
        plt.ylabel("PS Zero-point (mags)")
        pif CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_r_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_r_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_r_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_r_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

    if ("i" in zp_dict.keys()) or ("i" in zp_all_dict.keys()):

        # airmass correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "i" in zp_all_dict.keys():
            plt.errorbar(airmass_all_i,zp_all_i,yerr=zp_err_all_i,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "i" in zp_dict.keys():
            plt.errorbar(airmass_i,zp_i,yerr=zp_err_i,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Airmass (multiple of zenithal airmass)")
        plt.ylabel("PS Zero-point (mags)")
        plt.xlim(1.1,1.35)
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_airmass_i_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_airmass_i_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_airmass_i_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_airmass_i_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # timeobs correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "i" in zp_all_dict.keys():
            plt.errorbar(timeobs_all_i,zp_all_i,yerr=zp_err_all_i,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "i" in zp_dict.keys():
            plt.errorbar(timeobs_i,zp_i,yerr=zp_err_i,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Time of Observtion (hh.hh)")
        plt.ylabel("PS Zero-point (mags)")
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_i_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_i_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_i_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_i_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # FWHM correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "i" in zp_all_dict.keys():
            plt.errorbar(FWHM_all_i,zp_all_i,yerr=zp_err_all_i,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "i" in zp_dict.keys():
            plt.errorbar(FWHM_i,zp_i,yerr=zp_err_i,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Mean PSF FWHM (pixels)")
        plt.ylabel("PS Zero-point (mags)")
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_i_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_i_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_i_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_i_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

    if ("z" in zp_dict.keys()) or ("z" in zp_all_dict.keys()):

        # airmass correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "z" in zp_all_dict.keys():
            plt.errorbar(airmass_all_z,zp_all_z,yerr=zp_err_all_z,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "z" in zp_dict.keys():
            plt.errorbar(airmass_z,zp_z,yerr=zp_err_z,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Airmass (multiple of zenithal airmass)")
        plt.ylabel("PS Zero-point (mags)")
        plt.xlim(1.1,1.35)
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_airmass_z_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_airmass_z_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_airmass_z_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_airmass_z_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # timeobs correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "z" in zp_all_dict.keys():
            plt.errorbar(timeobs_all_z,zp_all_z,yerr=zp_err_all_z,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "z" in zp_dict.keys():
            plt.errorbar(timeobs_z,zp_z,yerr=zp_err_z,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Time of Observation (hh.hh)")
        plt.ylabel("PS Zero-point (mags)")
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_z_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_z_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_z_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_timeobs_z_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

        # FWHM correlation
        fig = plt.figure(1,figsize=(11,8.5))
        ax = fig.add_subplot(111)
        # all
        if "z" in zp_all_dict.keys():
            plt.errorbar(FWHM_all_z,zp_all_z,yerr=zp_err_all_z,fmt="o",color="blue",markersize=8,alpha=0.5)
        # cleaned
        if "z" in zp_dict.keys():
            plt.errorbar(FWHM_z,zp_z,yerr=zp_err_z,fmt="o",color="red",markersize=8,alpha=0.5)
        plt.xlabel("Mean PSF FWHM (pixels)")
        plt.ylabel("PS Zero-point (mags)")
        if CHECKZP==False:
            plt.ylim(22.,27.)
        elif CHECKZP==True:
            plt.ylim(-5.,5.)
        if exp=="LONG":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_z_longexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_z_longexp_calibrated.pdf",bbox_inches="tight")
        elif exp=="SHORT":
            if CHECKZP==False:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_z_shortexp.pdf",bbox_inches="tight")
            elif CHECKZP==True:
                plt.savefig(plotdir+"PS_zp_vs_FWHM_z_shortexp_calibrated.pdf",bbox_inches="tight")
        plt.close(fig)

def round_to_1(x):
    '''
    Rounds a float to one significant figure.
    '''
    if abs(x)>0.0:
        x_rounded = round(x, -int(floor(log10(abs(x)))))
    else:
        x_rounded = 0.0
    
    return x_rounded

def round_with_uncertainties(x,xerr):
    '''
    Truncate a float to n significant figures. May produce overflow in very last decimal place when q < 1.
    This can be removed by an extra formatted print.
    
    Input
    x : a float
    n : desired number of significant figures
    
    Output
    Float with only n s.f. and trailing zeros, but with a possible small overflow.
    '''
    
    if abs(xerr)>0.0:
        sgn=np.sign(x)
        x=abs(x)
        xerr_rounded = round_to_1(xerr) # overwrite input x
        n=int(np.log10(xerr_rounded/10.))
        x_rounded = sgn*round(x,abs(n))
    else:
        n = 2
        x_rounded = round(x,abs(n))
    
    return x_rounded, xerr

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
   return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def fgoodzp(zeropoint,passband,exptime,CHECKZP):
    '''
    Checks to see if the magnitude zeropoint obeys the lower-limit defined in the overlapconfig.yml configuration file.
    '''

    # regular calibration
    if CHECKZP==False:
        if passband=="i":
            zpmin = zpmin_i
            zpmax = zpmax_i
        elif passband=="r":
            if exptime==5.:
                zpmin = zpmin_r_short
                zpmax = zpmax_r_short
            else:
                zpmin = zpmin_r
                zpmax = zpmax_r
        elif passband=="z":
            zpmin = zpmin_z
            zpmax = zpmax_z
            
        if (zeropoint>=zpmin) and (zeropoint<=zpmax) and (zeropoint!=None):
            GOODZP=True
        else:
            GOODZP=False
    
    # check calibration
    if CHECKZP==True:
        if (zeropoint>=-1.) and (zeropoint<=1.0) and (zeropoint!=None):
            GOODZP=True
        else:
            GOODZP=False
    
    return GOODZP

def fpanstarrshist(panstarrsdata,plotdir):
    '''
    Creates a histogram of Pan-STARRS magnitudes.
    
    
    panstarrsdata[0] = r magnitude (PSF)
    panstarrsdata[1] = error in r magnitude (PSF)
    panstarrsdata[2] = r magnitude (Kron)
    panstarrsdata[3] = error in r magnitude (Kron)
    panstarrsdata[4] = r magnitude (Aperture)
    panstarrsdata[5] = error in r magnitude (Aperture)
    panstarrsdata[6] = i magnitude (PSF)
    panstarrsdata[7] = error in i magnitude (PSF)
    panstarrsdata[8] = i magnitude (Kron)
    panstarrsdata[9] = error in i magnitude (Kron)
    panstarrsdata[10] = i magnitude (Aperture)
    panstarrsdata[11] = error in i magnitude (Aperture)
    panstarrsdata[12] = z magnitude (PSF)
    panstarrsdata[13] = error in z magnitude (PSF)
    panstarrsdata[14] = z magnitude (Kron)
    panstarrsdata[15] = error in z magnitude (Kron)
    panstarrsdata[16] = z magnitude (Aperture)
    panstarrsdata[17] = error in z magnitude (Aperture)
    '''    
    
    # error-magnitude plots
    
    ########## PSF ##########
    
    print "plotting PSF magerr vs mag: r"
    # PSF magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[0],panstarrsdata[1],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_PSF_r.pdf",bbox_inches="tight")
    plt.close(fig)
    
    print "plotting PSF magerr vs mag: i"
    # PSF magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[6],panstarrsdata[7],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_PSF_i.pdf",bbox_inches="tight")
    plt.close(fig)
    
    print "plotting PSF magerr vs mag: z"
    # PSF magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[12],panstarrsdata[13],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_PSF_z.pdf",bbox_inches="tight")
    plt.close(fig)
    
    ########## KRON ##########
    
    print "plotting Kron magerr vs mag: r"
    # Kron magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[2],panstarrsdata[3],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_Kron_r.pdf",bbox_inches="tight")
    plt.close(fig)
    
    print "plotting Kron magerr vs mag: i"
    # Kron magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[8],panstarrsdata[9],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_Kron_i.pdf",bbox_inches="tight")
    plt.close(fig)
    
    print "plotting Kron magerr vs mag: z"
    # Kron magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[14],panstarrsdata[15],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_Kron_z.pdf",bbox_inches="tight")
    plt.close(fig)
    
    ########## APERTURE ##########
    
    print "plotting Aperture magerr vs mag: r"
    # Aperture magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[4],panstarrsdata[5],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_Aperture_r.pdf",bbox_inches="tight")
    plt.close(fig)
    
    print "plotting Aperture magerr vs mag: i"
    # Aperture magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[10],panstarrsdata[11],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_Aperture_i.pdf",bbox_inches="tight")
    plt.close(fig)
    
    print "plotting Aperture magerr vs mag: z"
    # Aperture magerr vs mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.hist2d(panstarrsdata[16],panstarrsdata[17],bins=100,cmin=1)
    plt.xlabel("PS Magnitude")
    plt.ylabel("PS Magnitude Uncertainty")
    plt.savefig(plotdir+"PS_allmagerr_Aperture_z.pdf",bbox_inches="tight")
    plt.close(fig)

def ftimeobs_to24h(time_fits):
    '''
    Converts the TIME-OBS FITS header keyword from hh:mm:ss.ss to hh.hh (i.e., to 24-hour clock time).
    '''

    time_fits_split = time_fits.split(":")
    time_24h        = np.float(time_fits_split[0]) + np.float(time_fits_split[1])/60. + np.float(time_fits_split[2])/3600.

    return time_24h

def fcheckobstime(passband,folder,timeobs,exptime):
    '''
    Checks to see if the observing time is satisfied for a given passband on a particular observing night.
    '''

    CHECKOBSTIME = True

    # 20111221
    if folder=="20111221":
        #   120s r                              120s i                                                                     120s z
        if ((timeobs>2.5) and (timeobs<3.0)) or ((timeobs>2.0) and (timeobs<2.5)) or ((timeobs>3.75) and (timeobs<4.5)) or ((timeobs>1.5) and (timeobs<2.0)) or ((timeobs>3.0) and (timeobs<4.0)):
            CHECKOBSTIME=False

    # 20111226
    elif folder=="20111226":
        #   120s r 
        if ((timeobs>2.0) and (timeobs<2.75)):
            CHECKOBSTIME=False

    # 20111227
    elif folder=="20111227":
        #  120s r                                                                      120s i                                                                    # 120s z
        if ((timeobs>2.0) and (timeobs<2.75)) or ((timeobs>5.5) and (timeobs<6.0)) or ((timeobs>3.25) and (timeobs<4.0)) or ((timeobs>5.0) and (timeobs<5.5)) or ((timeobs>2.75) and (timeobs<3.5)) or ((timeobs>4.5) and (timeobs<5.0)):
            CHECKOBSTIME=False

    # 20111228
    elif folder=="20111228":
        #  120s r                                                                                                           120s i                                # 120s z
        if ((timeobs>2.0) and (timeobs<2.75)) or ((timeobs>3.75) and (timeobs<4.5)) or ((timeobs>5.5) and (timeobs<6.0)) or ((timeobs>1.5) and (timeobs<2.25)) or ((timeobs>2.75) and (timeobs<3.5)):
            CHECKOBSTIME=False

    # 20120119
    elif folder=="20120119":
        #  120s r                                                                      120s i                                                                                                          # 120s z
        if ((timeobs>3.75) and (timeobs<4.5)) or ((timeobs>5.75) and (timeobs<7.5)) or ((timeobs>2.0) and (timeobs<2.75)) or ((timeobs>3.5) and (timeobs<4.0)) or ((timeobs>6.5) and (timeobs<7.0)) or ((timeobs>3.0) and (timeobs<3.75)) or ((timeobs>6.0) and (timeobs<7.0)):
            CHECKOBSTIME=False

    # 20120120
    elif folder=="20120120":
        #  120s r                                                                     # 120s z
        if ((timeobs>5.75) and (timeobs<6.5)) or ((timeobs>4.0) and (timeobs<4.5)) or ((timeobs>3.0) and (timeobs<4.0)):
            CHECKOBSTIME=False

    # 20120122
    elif folder=="20120122":
        CHECKOBSTIME=False

    # 20120123
    elif folder=="20120123":
        CHECKOBSTIME=False

    # 20120129
    elif folder=="20120129":
        # 120s z
        if ((timeobs>1.5) and (timeobs<2.0)):
            CHECKOBSTIME=False

    return CHECKOBSTIME
    

