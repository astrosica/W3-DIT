#!/usr/local/bin/python

#########################################################################################################
#                                                                                                       #
#                                 IMPORT CONFIGURATION FILE                                             #
#                                                                                                       #
#########################################################################################################


import yaml

with open("config.yml","r") as fconfig:
    config_data = yaml.safe_load(fconfig)
with open("overlapconfig.yml","r") as fconfig:
    overlapconfig_data = yaml.safe_load(fconfig)

parentdir         = config_data["parentdir"]
plotdir           = parentdir+"plots/"
passbands         = config_data["passbands"].split(" ")                # list of passbands to be used
calplotdir        = plotdir+"calibration/"                             # directory where calibration plots are stored
finalplotdir      = plotdir+"calibrated/"                              # directory where calibrated plots are stored
sexplotdir_cal    = calplotdir+"sextractor/"                           # subdirectory where SExtractor calibration plots are stored
apassplotdir_cal  = calplotdir+"apass/"                                # subdirectory where APASS calibration plots are stored
phottype          = overlapconfig_data["phottype"]                     # SExtractor photometry type; one of ISO, ISOCOR, or AUTO
detect_thresh     = float(overlapconfig_data["detect_thresh"])         # single pixel SNR threshold above RMS background
posmatch          = float(overlapconfig_data["posmatch"])/3600.        # maximum pointing offset for catalogue matching                [deg]
posmatch_lower    = float(overlapconfig_data["posmatch_lower"])/3600.  # lower-limit pointing offset for catalogue matching            [deg]
posmatch_upper    = float(overlapconfig_data["posmatch_upper"])/3600.  # upper-limit pointing offset for catalogue matching            [deg]
magmin            = float(overlapconfig_data["magmin"])                # minimum APASS magnitude                                       [deg]
magmax            = float(overlapconfig_data["magmax"])                # maximum APASS magnitude                                       [deg]
flagmax           = float(overlapconfig_data["flagmax"])               # maximum Source Extractor internal flag
clip_sigma        = float(overlapconfig_data["clip_sigma"])            # # of standard deviations for upper and lower clipping limits
zpmin_i           = float(overlapconfig_data["zpmin_i"])               # minimum limit on magnitude zeropoint in the i-band images
zpmin_r           = float(overlapconfig_data["zpmin_r"])               # minimum limit on magnitude zeropoint in the r-band images
sexcatdir         = config_data["sexcatdir"]                           # directory where SExtractor catalogue files are stored

#########################################################################################################
#                                                                                                       #
#                                         IMPORT PACKAGES                                               #
#                                                                                                       #
#########################################################################################################


import os
import numpy as np
from scipy import spatial
from astropy.io import ascii
from math import log10, floor
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit
from astropy.io.fits import getheader

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
from filefuncts import findfolderdir, appenditemtodict

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
    params_name = "default.param"
    nnw_name = "default.nnw"
    conv_name = "default.conv"
    catalog_name = "default.cat"

    if VERBOSE:
        verbose_type = "NORMAL"
    else:
        verbose_type = "QUIET"
    
    fp = open(dir+sextractor_config_name, "w")
    fp.write(sextractor_config.format(detect_thresh=detect_thresh, filter_name=conv_name,
        parameters_name=params_name, starnnw_name=nnw_name, verbose_type=verbose_type))
    fp.close()
    
    fp = open(dir+params_name, "w")
    fp.write(sextractor_params)
    fp.close()
    if VERBOSE:
        print "wrote file "+dir+params_name
    
    fp = open(dir+conv_name, "w")
    fp.write(default_conv)
    fp.close()
    if VERBOSE:
        print "wrote file "+dir+conv_name
    
    fp = open(dir+nnw_name, "w")
    fp.write(default_nnw)
    fp.close()
    if VERBOSE:
        print "wrote file "+dir+nnw_name
    
    return sextractor_config_name,params_name,nnw_name,conv_name,catalog_name

def sexcall(f,folder,fpath,objdir,catdir,bgdir,sexbool):
    '''
    Run source extractor on f, producing a catalogue as well as object and background maps
    f	    : file name to run source extractor on
    fpath   : directory path to f
    objdir  : directory to save object map in
    catdir  : directory to save catalogue in
    bgdir   :  directory to save background map in
    Returns location of output catalog
    '''

    # Split to make file name for catalogue, object map and background map filenames
    fname = f.split(".fts")[0]
    if sexbool==True:
        # Construct source extractor calls
        objsexcall = "sex -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME "+objdir+folder+"/"+fname+"_objects.fts "+fpath+f
        baksexcall = "sex -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME "+bgdir+folder+"/"+fname+"_background.fts "+fpath+f
        os.system(objsexcall)
        os.system(baksexcall)
    
    return catdir+folder+"/"+fname+".cat"

def fmagdiff(mag_sex,magerr_sex,mag_apass,magerr_apass):
    '''
    Calculates (APASS mag - SExtractor mag) and associated uncertainties.
    '''
    
    magdiff        = mag_apass - mag_sex
    magdiff_error  = np.sqrt((magerr_apass)**2. + (magerr_sex)**2.)
    
    return magdiff,magdiff_error

def magcalibration(file,folder,folderdir,apassdata,sexcat,overlapannfolderdir,VERBOSE,GENERATE):
    '''
    Builds KDTree of APASS catalog to be compared with Source Extractor catalogue.
    '''
    
    # load SExtractor catalog data
    sexdata      = ascii.read(sexcat)
    x_world      = sexdata["X_WORLD"]
    errx2_world  = sexdata["ERRX2_WORLD"]
    y_world      = sexdata["Y_WORLD"]
    erry2_world  = sexdata["ERRY2_WORLD"]
    
    # extract SExtractor R.A. and Dec.
    sexra        = x_world
    sexdec       = y_world
    sexradec     = [sexra,sexdec]
    sexradec     = np.array(np.transpose(sexradec))
    # extract APASS R.A. and Dec.
    apassra      = apassdata[:,0]
    apassdec     = apassdata[:,2]
    apassradec   = [apassra,apassdec]
    apassradec   = np.array(np.transpose(apassradec))
    
    # create APASS KDTree and use it to query SExtractor sources
    apasstree    = spatial.KDTree(apassradec)  # create KD-Tree
    matches      = apasstree.query(sexradec)   # find which APASS sources are closest to each DIT source
    dist         = np.array(matches[0])        # distance between nearest neighbours in degrees
    indices      = np.array(matches[1])        # APASS indices that match with SExtractor sources
    
    if VERBOSE:
        print len(apassradec), "stars extracted from APASS catalogue"
        print len(sexradec),   "stars found using Source Extractor"
        print len(indices),    "matching stars"
    
    if GENERATE:
        sexplotfolderdir    = sexplotdir_cal+folder+"/"
        apassplotfolderdir  = apassplotdir_cal+folder+"/"
        # histogram of pointing offsets between APASS and SExtractor catalogs
        pointingoffsethist(file,apassra,apassdec,sexra,sexdec,indices,dist,sexplotfolderdir)

    # clean up APASS overlap
    zeroslope,zeropoint,zeropoint_error,checkzp = fphotmatches(file,folder,folderdir,apassdata,indices,dist,sexcat,sexplotfolderdir,apassplotfolderdir,overlapannfolderdir,VERBOSE,GENERATE)
    
    return zeroslope,zeropoint,zeropoint_error,checkzp

def fphotmatches(
    file,folder,folderdir,apassdata,indices,dist,sexcat,sexplotfolderdir,apassplotfolderdir,overlapannfolderdir,VERBOSE,GENERATE,magmin=magmin,magmax=magmax,flagmax=flagmax):
    '''
    Cleans up APASS and Source Extractor photometry data using information stored in the APASS overlap.yml configuration file.
    '''

    # SExtractor catalog data
    sexdata           = ascii.read(sexcat)
    alpha_j2000       = sexdata["ALPHA_J2000"]           # Right Ascension (RA) of barycenter               [hh:mm:ss]
    x_world           = sexdata["X_WORLD"]               # barycenter position along world x-axis           [deg]
    errx2_world       = sexdata["ERRX2_WORLD"]           # variance of position along world x-axis          [deg**2]
    xerr_world        = np.sqrt(errx2_world)             # RMS uncertainty of position along world x-axis   [deg]
    delta_j2000       = sexdata["DELTA_J2000"]           # Declination (Dec) of barycenter                  [dd:mm:ss]
    y_world           = sexdata["Y_WORLD"]               # barycenter position along world y-axis           [deg]
    erry2_world       = sexdata["ERRY2_WORLD"]           # variance of position along world y-axis          [deg**2]
    yerr_world        = np.sqrt(erry2_world)             # RMS uncertainty of position along world y-axis   [deg]
    mag_iso           = sexdata["MAG_ISO"]               # isophotal manitude                               [mag]
    magerr_iso        = sexdata["MAGERR_ISO"]            # RMS uncertainty for ISO magnitude                [mag]
    flux_iso          = sexdata["FLUX_ISO"]              # flux density of ISO magnitude                    [ADU]
    fluxerr_iso       = sexdata["FLUX_ISO"]              # RMS uncertainty of ISO flux density              [ADU]
    mag_isocor        = sexdata["MAG_ISOCOR"]            # corrected isophotal magnitude                    [mag]
    magerr_isocor     = sexdata["MAGERR_ISOCOR"]         # RMS uncertainty for ISOCOR magnitude             [mag]
    flux_isocor       = sexdata["FLUX_ISOCOR"]           # flux density of ISOCOR magnitude                 [ADU]
    fluxerr_isocor    = sexdata["FLUX_ISOCOR"]           # RMS uncertainty of ISOCOR flux density           [ADU]
    mag_auto          = sexdata["MAG_AUTO"]              # kron-like elliptical aperture magnitude          [mag]
    magerr_auto       = sexdata["MAGERR_AUTO"]           # RMS uncertainty for AUTO magnitude               [mag]
    flux_auto         = sexdata["FLUX_AUTO"]             # flux density of AUTO magnitude                   [ADU]
    fluxerr_auto      = sexdata["FLUX_AUTO"]             # RMS uncertainty of AUTO flux density             [ADU]
    flags             = sexdata["FLAGS"]                 # extraction flags
    class_star        = sexdata["CLASS_STAR"]            # S/G classification
    
    header            = getheader(folderdir+file)
    passband          = header["FILTER"]
     
    # APASS catalogue data
    ra                = apassdata[:,0]                   # right ascension                                  [deg]
    ra_err            = apassdata[:,1]                   # error in right ascension                         [arcsec]
    dec               = apassdata[:,2]                   # declination                                      [deg]
    dec_err           = apassdata[:,3]                   # error in declination                             [arcsec]
    Nobs              = apassdata[:,4]                   # number of observations
    V_mag             = apassdata[:,5]                   # Johnson V magnitude                              [mag]
    V_magerr          = apassdata[:,6]                   # error in Johnson V magnitude                     [mag]
    B_mag             = apassdata[:,7]                   # Johnson B magnitude                              [mag]
    B_magerr          = apassdata[:,8]                   # error in Johnson B magnitude                     [mag]
    g_mag             = apassdata[:,9]                   # Sloan g' magnitude                               [mag]
    g_magerr          = apassdata[:,10]                  # error in Sloan g' magnitude                      [mag]
    r_mag             = apassdata[:,11]                  # Sloan r' magnitude                               [mag]
    r_magerr          = apassdata[:,12]                  # error in Sloan r' magnitude                      [mag]
    i_mag             = apassdata[:,13]                  # Sloan i' magnitude                               [mag]
    i_magerr          = apassdata[:,14]                  # error in Sloan i' magnitude                      [mag]
    
    # APASS matches with SExtractor sources
    ra_matches        = ra[indices]                      # right ascension                                  [deg]
    ra_err_matches    = ra_err[indices]                  # error in right ascension                         [arcsec]
    dec_matches       = dec[indices]                     # declination                                      [deg]
    dec_err_matches   = dec_err[indices]                 # error in declination                             [arcsec]
    Nobs_matches      = Nobs[indices]                    # number of observations
    V_mag_matches     = V_mag[indices]                   # Johnson V magnitude                              [deg]
    V_magerr_matches  = V_magerr[indices]                # error in Johnson V magnitude                     [deg]
    B_mag_matches     = B_mag[indices]                   # Johnson B magnitude                              [deg]
    B_magerr_matches  = B_magerr[indices]                # error in Johnson B magnitude                     [deg]
    g_mag_matches     = g_mag[indices]                   # Sloan g' magnitude                               [deg]
    g_magerr_matches  = g_magerr[indices]                # error in Sloan g' magnitude                      [deg]
    r_mag_matches     = r_mag[indices]                   # Sloan r' magnitude                               [deg]
    r_magerr_matches  = r_magerr[indices]                # error in Sloan r' magnitude                      [deg]
    i_mag_matches     = i_mag[indices]                   # Sloan i' magnitude                               [deg]
    i_magerr_matches  = i_magerr[indices]                # error in Sloan i' magnitude                      [deg]
    
    # obtain passband information
    if passband=="i":
        mag_apass     = i_mag_matches
        magerr_apass  = i_magerr_matches
    elif passband=="r":
        mag_apass     = r_mag_matches
        magerr_apass  = r_magerr_matches

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
    
    # clean up photometry catalog
    conditions = np.array(
        (dist<posmatch) & (flux>0.0) & (fluxerr>0.0) & (mag_apass<=magmax) & (mag_apass>=magmin) & (flags<=flagmax)
                           )
    i = np.array(np.where(conditions)[0])
    
    dist_clean            = dist[i]
    indices_clean         = indices[i]
    
    # APASS
    ra_clean              = ra_matches[i]
    ra_err_clean          = ra_err_matches[i]
    dec_clean             = dec_matches[i]
    dec_err_clean         = dec_err_matches[i]
    Nobs_clean            = Nobs_matches[i]
    r_mag_clean           = r_mag_matches[i]
    r_magerr_clean        = r_magerr_matches[i]
    i_mag_clean           = i_mag_matches[i]
    i_magerr_clean        = i_magerr_matches[i]
    mag_apass_clean       = mag_apass[i]
    magerr_apass_clean    = magerr_apass[i]
    
    # SOURCE EXTRACTOR
    alpha_j2000_clean     = alpha_j2000[i]
    x_world_clean         = x_world[i]
    xerr_world_clean      = xerr_world[i]
    delta_j2000_clean     = delta_j2000[i]
    y_world_clean         = y_world[i]
    yerr_world_clean      = yerr_world[i]
    mag_iso_clean         = mag_iso[i]
    magerr_iso_clean      = magerr_iso[i]
    flux_iso_clean        = flux_iso[i]
    fluxerr_iso_clean     = fluxerr_iso[i]
    mag_isocor_clean      = mag_isocor[i]
    magerr_isocor_clean   = magerr_isocor[i]
    flux_isocor_clean     = flux_isocor[i]
    fluxerr_isocor_clean  = fluxerr_isocor[i]
    mag_auto_clean        = mag_auto[i]
    magerr_auto_clean     = magerr_auto[i]
    flux_auto_clean       = flux_auto[i]
    fluxerr_auto_clean    = fluxerr_auto[i]
    flags_clean           = flags[i]
    class_star_clean      = class_star[i]
    mag_sex_clean         = mag_sex[i]
    magerr_sex_clean      = magerr_sex[i]

    # perform sigma-clipping on magnitude zero points
    ra_sc,dec_sc,x_world_sc,y_world_sc,mag_apass_sc,magerr_apass_sc,mag_sex_sc,magerr_sex_sc,zp_sc=fzeropoint_sc(ra_clean,dec_clean,x_world_clean,y_world_clean,mag_apass_clean,magerr_apass_clean,mag_sex_clean,magerr_sex_clean,clip_sigma)
    
    # perform linear fitting to calculate magnitude zero point and slope
    zeroslope,zeropoint,zeropoint_error = fitzeropoint(file,mag_sex_clean,magerr_sex_clean,mag_apass_clean,magerr_apass_clean,sexplotfolderdir)
    zeroslope_sc = fitzeropoint(file,mag_sex_sc,magerr_sex_sc,mag_apass_sc,magerr_apass_sc,sexplotfolderdir,sc=True)
    
    checkzp = fcheckzp(zeropoint,passband)
    
    if GENERATE:
        # create annotation files for APASS matches
        '''
        writeannfiles(file,ra_clean,dec_clean,"_APASS_clean",overlapannfolderdir,VERBOSE)                # cleaned APASS pointings
        writeannfiles(file,x_world_clean,y_world_clean,"_Sextractor_clean",overlapannfolderdir,VERBOSE)  # cleaned SExtractor pointings
        writeannfiles(file,ra_sc,dec_sc,"_APASS_sc",overlapannfolderdir,VERBOSE)                         # cleaned and sigma-clipped
        writeannfiles(file,x_world_sc,y_world_sc,"_SExtractor_sc",overlapannfolderdir,VERBOSE)           # cleaned and sigma-clipped

        # plot ISO, ISOCOR, and AUTO magnitudes
        comparemags(file,mag_iso_clean,magerr_iso_clean,mag_isocor_clean,magerr_isocor_clean,mag_auto_clean,magerr_auto_clean,mag_apass_clean,magerr_apass_clean,sexplotfolderdir)

        # plot SExtractor and APASS positional errors versus SExtractor's ISO, ISOCOR, and AUTO magnitudes
        poserrvsmags(file,ra_err_clean,dec_err_clean,xerr_world_clean,yerr_world_clean,mag_iso_clean,magerr_iso_clean,mag_isocor_clean,magerr_isocor_clean,mag_auto_clean,magerr_auto_clean,mag_apass_clean,magerr_apass_clean,sexplotfolderdir)

        # plot (APASS mag - SExtractor mag) vs. APASS mag
        plotmagsoln(file,mag_apass_clean,magerr_apass_clean,mag_sex_clean,magerr_sex_clean,sexplotfolderdir)
        #plotmagsoln(file,mag_apass_sc,magerr_apass_sc,mag_sex_sc,magerr_sex_sc,sexplotfolderdir,sc=True)

        # plot histograms of (APASS mag - SExtractor mag)
        plotmagsolnhist(file,mag_apass_clean,mag_sex_clean,dist_clean,sexplotfolderdir)

        # plot histograms of APASS and SExtractor mag uncertainties
        plotmagerrorhist(file,magerr_apass_clean,magerr_iso_clean,magerr_isocor_clean,magerr_auto_clean,sexplotfolderdir)
        '''
        
        if (zeroslope==True) and (checkzp==True):
            # color-magnitude diagram using APASS catalogue only
            apasscmdiagram(file,ra_clean,dec_clean,r_mag_clean,r_magerr_clean,i_mag_clean,i_magerr_clean,apassplotfolderdir)
        
    return zeroslope,zeropoint,zeropoint_error,checkzp

def fcheckzp(zeropoint,passband):
    '''
    Checks to see if the magnitude zeropoint obeys the lower-limit defined in the overlapconfig.yml configuration file.
    '''
    
    if passband=="i":
        zpmin = zpmin_i
    elif passband=="r":
        zpmin = zpmin_r
        
    if zeropoint>=zpmin:
        checkzp=True
    else:
        checkzp=False
    
    return checkzp

def sexcmdiagram(file_i,folder_i,file_r,folder_r,catdir,zeropoints,plotdir,addbirs=True):
    '''
    Plots the colour-magnitude (CM) digram of r vs. (r - i) for calibrated SExtractor magnitudes.
    '''
    
    sexcatfile_i = catdir+folder_i+"/"+file_i.replace(".fts",".cat")
    sexcatfile_r = catdir+folder_r+"/"+file_r.replace(".fts",".cat")
    
    folderdir_i  = findfolderdir(folder_i)
    folderdir_r  = findfolderdir(folder_r)
    
    sexdata_i = ascii.read(sexcatfile_i)
    sexdata_r = ascii.read(sexcatfile_r)
    
    if phottype=="AUTO":
        mag_i,magerr_i = sexdata_i["MAG_AUTO"],sexdata_i["MAGERR_AUTO"]
        mag_r,magerr_r = sexdata_r["MAG_AUTO"],sexdata_r["MAGERR_AUTO"]
    elif phottype=="ISO":
        mag_i,magerr_i = sexdata_i["MAG_ISO"],sexdata_i["MAGERR_ISO"]
        mag_r,magerr_r = sexdata_r["MAG_ISO"],sexdata_r["MAGERR_ISO"]
    elif phottype=="ISOCOR":
        mag_i,magerr_i = sexdata_i["MAG_ISOCOR"],sexdata_i["MAGERR_ISOCOR"]
        mag_r,magerr_r = sexdata_r["MAG_ISOCOR"],sexdata_r["MAGERR_ISOCOR"]

    ra_i      = sexdata_i["X_WORLD"]
    dec_i     = sexdata_i["Y_WORLD"]
    flags_i   = sexdata_i["FLAGS"]
    radec_i   = [ra_i,dec_i]
    radec_i   = np.array(np.transpose(radec_i))
    
    ra_r      = sexdata_r["X_WORLD"]
    dec_r     = sexdata_r["Y_WORLD"]
    flags_r   = sexdata_r["FLAGS"]
    radec_r   = [ra_r,dec_r]
    radec_r   = np.array(np.transpose(radec_r))
    
    kdtree    = spatial.KDTree(radec_i)
    matches   = kdtree.query(radec_r)
    dist      = np.array(matches[0])
    indices   = np.array(matches[1])
    
    ra_i      = ra_i[indices]
    dec_i     = dec_i[indices]
    flags_i   = flags_i[indices]
    mag_i     = mag_i[indices]
    magerr_i  = magerr_i[indices]
    
    zp_zperr_i = zeropoints[folderdir_i+file_i][0]
    zp_zperr_r = zeropoints[folderdir_r+file_r][0]
    
    zp_i,zperr_i = zp_zperr_i[0],zp_zperr_i[1]
    zp_r,zperr_r = zp_zperr_r[0],zp_zperr_r[1]
    
    mag_i_cal,magerr_i_cal = mag_i+zp_i,np.sqrt((magerr_i)**2. + (zperr_i)**2.)
    mag_r_cal,magerr_r_cal = mag_r+zp_r,np.sqrt((magerr_r)**2. + (zperr_r)**2.)
    
    conditions = np.array( (dist<posmatch) & (mag_i_cal>0.0) & (mag_r_cal>0.0) & (flags_i<=flagmax) & (flags_r<=flagmax) )
    i = np.array(np.where(conditions)[0])
    
    mag_i_clean     = mag_i_cal[i]
    magerr_i_clean  = magerr_i_cal[i]
    mag_r_clean     = mag_r_cal[i]
    magerr_r_clean  = magerr_r_cal[i]
    
    ri_colour       = mag_r_clean - mag_i_clean
    ri_colour_error = np.sqrt((magerr_r_clean)**2. + (magerr_i_clean)**2.)
    
    if addbirs==True:
        '''
        Adds location of known BIRS.
        '''
        
        lines    = np.loadtxt("/mnt/raid-project/hp/campbell/DIT/annFiles/allbirs.txt",delimiter=" ")
        birs_ra  = lines[:,0]
        birs_dec = lines[:,1]
        
        # use KD-Tree to find matches between APASS sources and known BIRS
        sexradec   = [ra_r,dec_r]
        sexradec   = np.array(np.transpose(sexradec))
        birsradec  = [birs_ra,birs_dec]
        birsradec  = np.array(np.transpose(birsradec))
        birstree   = spatial.KDTree(birsradec)
        matches    = birstree.query(sexradec)
        dist       = np.array(matches[0])
        indices    = np.array(matches[1])
        
        conditions = np.array(dist<posmatch)
        
        i_birs = np.array(np.where(conditions)[0])

        #dist_birs           = dist[i_birs]
        indices_birs         = indices[i_birs]
        
        mag_i_birs            = mag_i_clean[indices_birs]
        magerr_i_birs         = magerr_i_clean[indices_birs]
        mag_r_birs            = mag_r_clean[indices_birs]
        magerr_r_birs         = magerr_r_clean[indices_birs]
        ri_colour_birs        = ri_colour[indices_birs]
        ri_colour_error_birs  = ri_colour_error[indices_birs]
        
        print len(indices_birs), "matches with known BIRS for APASS colour-magnitude diagram"
        
    if addbirs==True:
        newfname = "colourmag_r-i_BIRS.pdf"
    else:
        newfname = "colourmag_r-i.pdf"
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(ri_colour,mag_r_clean,xerr=ri_colour_error,yerr=magerr_r_clean,fmt="o",alpha=0.5,markersize=8,color="purple",label="SExtractor")
    if addbirs==True:
        plt.errorbar(ri_colour_birs,mag_r_birs,xerr=ri_colour_error_birs,yerr=magerr_r_birs,fmt="*",alpha=0.8,markersize=10,color="red",label="BIRS")
    #plt.legend()
    plt.gca().invert_yaxis()
    plt.xlabel("r-i Magnitude (mag)")
    plt.ylabel("SExtractor r Magnitude (mag)")
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)

def apasscmdiagram(file,ra,dec,r_mag,r_magerr,i_mag,i_magerr,plotdir,addbirs=True):
    '''
    Plots the colour-magnitude (CM) diagram of r vs. (r - i) for APASS magnitudes.
    '''
    
    conditions = np.array( (r_mag>0.0) & (i_mag>0.0) )

    i = np.array(np.where(conditions)[0])
    
    r_mag     = r_mag[i]
    r_magerr  = r_magerr[i]
    i_mag     = i_mag[i]
    i_magerr  = i_magerr[i]
    
    ri_colour        = r_mag - i_mag
    ri_colour_error  = np.sqrt((r_magerr)**2. + (i_magerr)**2.)
    
    if addbirs==True:
        '''
        Adds location of known BIRS.
        '''
        
        lines    = np.loadtxt("/mnt/raid-project/hp/campbell/DIT/annFiles/allbirs.txt",delimiter=" ")
        birs_ra  = lines[:,0]
        birs_dec = lines[:,1]
        
        # use KD-Tree to find matches between APASS sources and known BIRS
        sexradec   = [ra,dec]
        sexradec   = np.array(np.transpose(sexradec))
        birsradec  = [birs_ra,birs_dec]
        birsradec  = np.array(np.transpose(birsradec))
        birstree   = spatial.KDTree(birsradec)
        matches    = birstree.query(sexradec)
        dist       = np.array(matches[0])
        indices    = np.array(matches[1])
        
        conditions = np.array(dist<posmatch)
        
        i_birs = np.array(np.where(conditions)[0])

        #dist_birs           = dist[i_birs]
        indices_birs         = indices[i_birs]

        r_mag_birs            = r_mag[indices_birs]
        r_magerr_birs         = r_magerr[indices_birs]
        ri_colour_birs        = ri_colour[indices_birs]
        ri_colour_error_birs  = ri_colour_error[indices_birs]
        
        print len(indices_birs), "matches with known BIRS for APASS colour-magnitude diagram"
        
    if addbirs==True:
        newfname = file.replace(".fts","_r-i_BIRS.pdf")
    else:
        newfname = file.replace(".fts","_r-i.pdf")
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(ri_colour,r_mag,xerr=ri_colour_error,yerr=r_magerr,fmt="o",alpha=0.5,markersize=8,color="purple",label="APASS")
    if addbirs==True:
        plt.errorbar(ri_colour_birs,r_mag_birs,xerr=ri_colour_error_birs,yerr=r_magerr_birs,fmt="*",alpha=0.8,markersize=10,color="red",label="BIRS")
    #plt.legend()
    plt.gca().invert_yaxis()
    plt.xlabel("r-i Magnitude (mag)")
    plt.ylabel("APASS r Magnitude (mag)")
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)

def plotmagerrorhist(file,magerr_apass,magerr_iso,magerr_isocor,magerr_auto,plotdir):
    '''
    Plot histograms of APASS and SExtractor ISO, ISOCOR, and AUTO magnitude uncertainties.
    '''
    
    deltaerr = 0.05
    errlower = np.min([np.min(magerr_apass),np.min(magerr_iso),np.min(magerr_isocor),np.min(magerr_auto)])
    errupper = np.max([np.max(magerr_apass),np.max(magerr_iso),np.max(magerr_isocor),np.max(magerr_auto)])
    
    N = len(np.arange(errlower,errupper,deltaerr))
    commonparams = dict(bins=N,range=(errlower,errupper))
    histbins = np.linspace(errlower,errupper,N,endpoint=True)
    
    legend_apass = mpatches.Patch(color="blue",label="APASS")
    legend_iso = mpatches.Patch(color="green",label="ISO")
    legend_isocor = mpatches.Patch(color="red",label="ISOCOR")
    legend_auto = mpatches.Patch(color="turquoise",label="AUTO")
    
    newfname = file.replace(".fts","_magerr_hist.pdf")
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    pylab.hist((magerr_apass,magerr_iso,magerr_isocor,magerr_auto),**commonparams)
    plt.xlabel("Magnitude Uncertainty")
    plt.xlim(errlower,errupper)
    plt.ylabel("N")
    plt.legend(loc="upper right",handles=[legend_apass,legend_iso,legend_isocor,legend_auto],prop={"size":14})
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)

def fgauss(x,amp,mean,sigma):
    '''
    Gaussian distribution function.
    x      : input variable
    amp    : maximum amplitude
    mean   : distribution mean
    sigma  : standard deviation
    '''
    y = amp * np.exp(-((x - mean)**2. / (2. * sigma**2.)))
    return y

def round_to_1(x):
    '''
    Rounds a float to one significant figure.
    '''
    
    x_rounded = round(x, -int(floor(log10(abs(x)))))
    
    return x_rounded

def round_with_uncertainties(x,xerr):
    '''
    Truncate a float to n significant figures.  May produce overflow in 
    very last decimal place when q < 1.  This can be removed by an extra 
    formatted print. 
    Arguments:
      x : a float
      n : desired number of significant figures
    Returns:
    Float with only n s.f. and trailing zeros, but with a possible small overflow.
    '''
    
    sgn=np.sign(x)
    x=abs(x)
    xerr_rounded = round_to_1(xerr) # overwrite input x
    n=int(np.log10(xerr_rounded/10.))
    x_rounded = round(x,abs(n))
    
    return x_rounded,xerr_rounded
    

def plotmagsolnhist(file,mag_apass,mag_sex,dist,plotdir,deltamag=0.05,poslower=posmatch_lower,posupper=posmatch_upper):
    '''
    Plots histograms of (APASS mag - SExtractor mag) for matches within both x and y arcseconds defined in overlapconfig.yml.
    '''
    
    condition_poslower            = np.array(dist<=poslower)
    condition_posupper            = np.array( (dist>poslower) & (dist<=posupper) )
    i_poslower                    = np.array(np.where(condition_poslower)[0])
    i_posupper                    = np.array(np.where(condition_posupper)[0])
    mag_apass_poslower            = mag_apass[i_poslower]
    mag_apass_posupper            = mag_apass[i_posupper]
    mag_sex_poslower              = mag_sex[i_poslower]
    mag_sex_posupper              = mag_sex[i_posupper]
    magdiff_poslower              = mag_apass_poslower - mag_sex_poslower
    magdiff_posupper              = mag_apass_posupper - mag_sex_posupper
    maglower                      = np.min([np.min(magdiff_poslower),np.min(magdiff_poslower)])
    magupper                      = np.max([np.max(magdiff_posupper),np.max(magdiff_posupper)])
    histmin                       = np.min([np.min(magdiff_poslower),np.min(magdiff_posupper)])
    histmax                       = np.max([np.max(magdiff_poslower),np.max(magdiff_posupper)])
    N                             = len(np.arange(histmin,histmax,deltamag))
    commonparams                  = dict(bins=N,range=(histmin,histmax))
    histbins                      = np.linspace(histmin,histmax,N,endpoint=True)
    
    # fit Gaussian
    yhist_poslower,xhist_poslower  = np.histogram(magdiff_poslower,bins=histbins)
    yhist_posupper,xhist_posupper  = np.histogram(magdiff_posupper,bins=histbins)
    xhist_halfwidth_poslower       = 0.5*abs(xhist_poslower[1]-xhist_poslower[0])
    xhist_halfwidth_posupper       = 0.5*abs(xhist_posupper[1]-xhist_posupper[0])
    xhist_gauss_poslower           = xhist_poslower[:-1]+xhist_halfwidth_poslower
    xhist_gauss_posupper           = xhist_posupper[:-1]+xhist_halfwidth_posupper
    mu_poslower, sigma_poslower    = norm.fit(magdiff_poslower)
    mu_posupper, sigma_posupper    = norm.fit(magdiff_posupper)
    
    try:
        fit=True
        popt_poslower,pcov_poslower = curve_fit(
            fgauss,xhist_gauss_poslower,yhist_poslower,np.array([np.max(yhist_poslower),mu_poslower,sigma_poslower]),maxfev=1000)
        popt_poslower_error = np.sqrt(np.diag(pcov_poslower))
        mean_poslower = popt_poslower[1]
        mean_poslower_error = popt_poslower_error[1]
        mean_poslower_rounded,mean_poslower_error_rounded = round_with_uncertainties(mean_poslower,mean_poslower_error)
        popt_posupper,pcov_posupper = curve_fit(
            fgauss,xhist_gauss_posupper,yhist_posupper,np.array([np.max(yhist_posupper),mu_posupper,sigma_posupper]),maxfev=1000)
    except RuntimeError:
        fit=False
    
    newfname = file.replace(".fts","_magsoln_hist.pdf")
    
    legend_poslower = mpatches.Patch(color="blue",label=str(poslower*3600.)+" arcsec")
    legend_posupper = mpatches.Patch(color="green",label=str(posupper*3600.)+" arcsec")

    # histogram of (APASS mag - SExtractor mag)
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    # histograms
    pylab.hist((magdiff_poslower,magdiff_posupper),**commonparams)
    # fits
    if fit==True:
        plt.plot(np.linspace(histmin,histmax,100),fgauss(np.linspace(histmin,histmax,100),*popt_poslower),"b--",linewidth=2)
        plt.plot(np.linspace(histmin,histmax,100),fgauss(np.linspace(histmin,histmax,100),*popt_posupper),"g--",linewidth=2)
        plt.axvline(popt_poslower[1],color="blue",linestyle="dashed",linewidth=2)
        plt.axvline(popt_posupper[1],color="green",linestyle="dashed",linewidth=2)
        plt.text(0.05,0.7,"mean: "+str(mean_poslower_rounded)+" +/- "+str(mean_poslower_error_rounded),transform=ax.transAxes)
    elif fit==False:
        pass
    plt.xlabel("(APASS Magnitude - SExtractor Magnitude)")
    plt.xlim(histmin,histmax)
    plt.ylabel("N")
    plt.legend(loc="upper left",handles=[legend_poslower,legend_posupper],prop={"size":14})
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)
    
def poserrvsmags(file,ra_err,dec_err,xerr_world,yerr_world,mag_iso,magerr_iso,mag_isocor,magerr_isocor,mag_auto,magerr_auto,mag_apass,magerr_apass,plotdir):
    '''
    Plots Source Extractor ISO, ISOCOR, and AUTO positional uncertainty versus APASS magnitude.
    '''
    
    poserr_APASS   = np.sqrt((ra_err)**2. + (dec_err)**2.)
    poserr_SEX     = np.sqrt((xerr_world)**2. + (yerr_world)**2.)*3600.
    
    newfname_SEX   = file.replace(".fts","_poserrvsmags_sextractor.pdf")
    newfname_APASS = file.replace(".fts","_poserrvsmags_apass.pdf")
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(mag_iso,poserr_SEX,xerr=magerr_iso,fmt="o",color="purple",label="ISO")
    plt.errorbar(mag_isocor,poserr_SEX,xerr=magerr_isocor,fmt="o",color="orange",label="ISOCOR")
    plt.errorbar(mag_auto,poserr_SEX,xerr=magerr_auto,fmt="o",color="turquoise",label="AUTO")
    # plot log
    ax.set_yscale("log")
    plt.xlabel("Source Extractor Magnitude (mag)")
    plt.ylabel("Sextractor RMS Positional Uncertainty (arcseconds)")
    plt.legend(loc="upper left",numpoints=1)
    plt.savefig(plotdir+newfname_SEX,bbox_inches="tight")
    plt.close(fig)
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(mag_apass,poserr_APASS,xerr=magerr_apass,fmt="o",color="purple")
    # plot log
    ax.set_yscale("log")
    plt.xlabel("APASS Magnitude (mag)")
    plt.ylabel("APASS RMS Positional Uncertainty (arcseconds)")
    plt.savefig(plotdir+newfname_APASS,bbox_inches="tight")
    plt.close(fig)

def comparemags(file,mag_iso,magerr_iso,mag_isocor,magerr_isocor,mag_auto,magerr_auto,mag_apass,magerr_apass,plotdir):
    '''
    Plots Source Extractor's ISO, ISOCOR, and AUTO magnitudes versus APASS magnitudes.
    '''
    
    # APASS mag - SExtractor mag
    magdiff_iso,    magdiff_iso_error     = fmagdiff(mag_iso,magerr_iso,mag_apass,magerr_apass)
    magdiff_isocor, magdiff_isocor_error  = fmagdiff(mag_isocor,magerr_isocor,mag_apass,magerr_apass)
    magdiff_auto,   magdiff_auto_error    = fmagdiff(mag_auto,magerr_auto,mag_apass,magerr_apass)
    
    newfname = file.replace(".fts","_comparemags.pdf")
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    # (APASS mag - SExtractor mag) vs. APASS mag
    plt.errorbar(mag_apass,magdiff_iso,xerr=magerr_apass,yerr=magdiff_iso_error,fmt="o",color="purple",label="ISO")
    plt.errorbar(mag_apass,magdiff_isocor,xerr=magerr_apass,yerr=magdiff_isocor_error,fmt="o",color="orange",label="ISOCOR")
    plt.errorbar(mag_apass,magdiff_auto,xerr=magerr_apass,yerr=magdiff_auto_error,fmt="o",color="turquoise",label="AUTO")
    plt.xlabel("APASS Magnitude")
    plt.ylabel("Source Extractor Magnitude")
    plt.legend(loc="lower left",numpoints=1)
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)

def apasshist(apassdata,maglower,magupper,deltamag,plotdir):
    '''
    Creates a histogram of APASS magnitudes.
    
    apassdata[:,0]  = right ascension in degrees
    apassdata[:,1]  = error in right ascension in arcseconds
    apassdata[:,2]  = declination in degrees
    apassdata[:,3]  = error in declination in arcseconds
    apassdata[:,4]  = number of observations
    apassdata[:,5]  = Johnson V magnitude
    apassdata[:,6]  = error in Johnson V magnitude
    apassdata[:,7]  = Johnson B magnitude
    apassdata[:,8]  = error in Johnson B magnitude
    apassdata[:,9]  = Sloan g' magnitude
    apassdata[:,10] = error in Sloan g' magnitude
    apassdata[:,11] = Sloan r' magnitude
    apassdata[:,12] = error in Sloan r' magnitude
    apassdata[:,13] = Sloan i' magnitude
    apassdata[:,14] = error in Sloan i' magnitude
    '''
    
    magcols        = [5,7,9,11,13]
    magerrcols     = [6,8,10,12,14]
    magstrings     = ["Johnson V magnitude","Johnson B magnitude","Sloan g' magnitude","Sloan r' magnitude","Sloan i' magnitude"]
    hist_fnames    = ["Johnson_V_hist","Johnson_B_hist","Sloan_g_hist","Sloan_r_hist","Sloan_i_hist"]
    magerr_fnames  = ["Johnson_V_magerr","Johnson_B_magerr","Sloan_g_magerr","Sloan_r_magerr","Sloan_i_magerr"]
    colours        = ["purple","orange","green","yellow","turquoise"]
    
    # combined error-magnitude plot
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    
    for i in range(len(magcols)):
        magcol        = magcols[i]
        magerrcol     = magerrcols[i]
        mag_apass     = apassdata[:,magcol]
        magerr_apass  = apassdata[:,magerrcol]
        magstring     = magstrings[i]
        magerr_fname  = magerr_fnames[i]
        colour        = colours[i]
        
        plt.plot(mag_apass,magerr_apass,"ro",color=colour,alpha=0.5,label=magstring)
    
    plt.axhline(0,linestyle="dashed",color="black")
    plt.xlabel("Magnitude")
    plt.ylabel("Magnitude Uncertainty")
    plt.legend(loc="upper left",numpoints=1)
    plt.savefig(plotdir+"allmagerr.pdf",bbox_inches="tight")
    plt.close(fig)
    
    # combined magnitude histogram plot using input limits
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    
    hist1 = apassdata[:,magcols[0]]
    hist2 = apassdata[:,magcols[1]]
    hist3 = apassdata[:,magcols[2]]
    hist4 = apassdata[:,magcols[3]]
    hist5 = apassdata[:,magcols[4]]
    
    N = len(np.arange(maglower,magupper,deltamag))
    commonparams = dict(bins=N,range=(maglower,magupper))
    
    pylab.hist((hist1,hist2,hist3,hist4,hist5),**commonparams)
    
    hist1_legend = mpatches.Patch(color="blue",label=magstrings[0])
    hist2_legend = mpatches.Patch(color="green",label=magstrings[1])
    hist3_legend = mpatches.Patch(color="red",label=magstrings[2])
    hist4_legend = mpatches.Patch(color="turquoise",label=magstrings[3])
    hist5_legend = mpatches.Patch(color="purple",label=magstrings[4])

    plt.xlabel("Magnitude")
    plt.ylabel("N")
    plt.legend(loc="upper left",handles=[hist1_legend,hist2_legend,hist3_legend,hist4_legend,hist5_legend])
    plt.savefig(plotdir+"allhist_cut.pdf",bbox_inches="tight")
    plt.close(fig)
    
    # full magnitude histogram
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)

    hist1 = apassdata[:,magcols[0]]
    hist2 = apassdata[:,magcols[1]]
    hist3 = apassdata[:,magcols[2]]
    hist4 = apassdata[:,magcols[3]]
    hist5 = apassdata[:,magcols[4]]

    maglower = np.min([np.min(hist1),np.min(hist2),np.min(hist3),np.min(hist4),np.min(hist5)])
    magupper = np.max([np.max(hist1),np.max(hist2),np.max(hist3),np.max(hist4),np.max(hist5)])
    
    N = len(np.arange(maglower,magupper,deltamag))
    commonparams = dict(bins=N,range=(maglower,magupper))
    
    pylab.hist((hist1,hist2,hist3,hist4,hist5),**commonparams)
    
    hist1_legend = mpatches.Patch(color="blue",label=magstrings[0])
    hist2_legend = mpatches.Patch(color="green",label=magstrings[1])
    hist3_legend = mpatches.Patch(color="red",label=magstrings[2])
    hist4_legend = mpatches.Patch(color="turquoise",label=magstrings[3])
    hist5_legend = mpatches.Patch(color="purple",label=magstrings[4])

    plt.xlabel("Magnitude")
    plt.ylabel("N")
    plt.legend(loc="upper right",handles=[hist1_legend,hist2_legend,hist3_legend,hist4_legend,hist5_legend])
    plt.savefig(plotdir+"allhist.pdf",bbox_inches="tight")
    plt.close(fig)

def writeannfiles(file,ra_all,dec_all,fname_suffix,overlapannfolderdir,VERBOSE,colour="RED"):
    '''
    Writes APASS matches to annotation files in both kvis and ds9 formats.
    '''
    
    color_kvis = colour
    color_ds9 = color_kvis.lower()
    fname_kvis = file.replace(".fts",fname_suffix+".ann")
    fname_ds9 = file.replace(".fts",fname_suffix+".reg")
    f_kvis = open(overlapannfolderdir+fname_kvis,"w+")
    f_ds9 = open(overlapannfolderdir+fname_ds9,"w+")
    f_kvis.write("COLOR %s\n\n" % color_kvis)
    for i in range(len(ra_all)):
        ra_i = ra_all[i]
        dec_i = dec_all[i]
        pos_err = 0.002
        f_kvis.write("CIRCLE {0} {1} {2}\n\n".format(ra_i,dec_i,pos_err))
        f_ds9.write("fk5;circle {0} {1} {2} # color={3}\n\n".format(ra_i,dec_i,pos_err,color_ds9))
    f_kvis.close()
    f_ds9.close()

    if VERBOSE:
        print "wrote kvis annotation file "+fname_kvis
        print "wrote ds9 annotation file "+fname_ds9

def fitzeropoint(file,mag_sex,magerr_sex,mag_apass,magerr_apass,plotdir,sc=False,yerrors=False,plotmedian=True,plotweightedavg=False,cleanmags=True,magmin=11.,magmax=14.):
    '''
    Fits a straight line to (SExtractor - APASS) vs. APASS magnitudes
    to measure magnitude zero point and slope in the data.
    '''
    
    def fline(mag_apass,slope,zeropoint):
        y =  slope*mag_apass + zeropoint
        return y
    
    def fline_odr(params,mag_apass):
        slope = params[0]
        zeropoint = params[1]
        y = slope*mag_apass + zeropoint
        return y
    
    conditions = np.array( (mag_apass>=magmin) & (mag_apass<=magmax) )
    i = np.array(np.where(conditions)[0])
    
    mag_sex_clean       = mag_sex[i]
    magerr_sex_clean    = magerr_sex[i]
    mag_apass_clean     = mag_apass[i]
    magerr_apass_clean  = magerr_apass[i]
    
    if cleanmags==True:
        mag_sex_raw,magerr_sex_raw,mag_apass_raw,magerr_apass_raw = mag_sex.copy(),magerr_sex.copy(),mag_apass.copy(),magerr_apass.copy()
        mag_sex,magerr_sex,mag_apass,magerr_apass                 = mag_sex_clean.copy(),magerr_sex_clean.copy(),mag_apass_clean.copy(),magerr_apass_clean.copy()
    
    magdiff,magdiff_error = fmagdiff(mag_sex,magerr_sex,mag_apass,magerr_apass)
    magdiff_median        = np.median(magdiff)
    magdiff_median_error  = np.std(magdiff)/np.sqrt(len(magdiff))
    magdiff_weights       = 1./(magdiff_error)**2.
    magdiff_weightedavg   = np.average(magdiff,weights=magdiff_weights)
    
    try:
        curvefit=True
        guess = np.array([0.0,25.5])
        if yerrors==False:
            popt,pcov           = curve_fit(fline,mag_apass,magdiff,p0=guess,maxfev=1000)
            popt_uncertainties  = np.sqrt(np.diag(pcov))
        elif yerrors==True:
            popt,pcov            = curve_fit(fline,mag_apass,magdiff,p0=guess,sigma=magdiff_error,maxfev=1000)
            popt_uncertainties   = np.sqrt(np.diag(pcov))
        maglower                 = np.min(mag_apass)
        magupper                 = np.max(mag_apass)
        mag_apass_fit            = np.linspace(maglower,magupper,100)
        deltamag_fit             = fline(mag_apass_fit,*popt)
        
        slope                    = popt[0]
        slope_error              = popt_uncertainties[0]
        zeropoint                = popt[1]
        zeropoint_error          = popt_uncertainties[1]
        
        slope_rounded,slope_error_rounded = round_with_uncertainties(slope,slope_error)
        zeropoint_rounded,zeropoint_error_rounded = round_with_uncertainties(zeropoint,zeropoint_error)
        magdiff_median_rounded,magdiff_median_error_rounded = round_with_uncertainties(magdiff_median,magdiff_median_error)
        
        errorsig = 3.0
        if slope<0.0:
            if (slope+errorsig*slope_error)>=0.0:
                zeroslope=True
            else:
                zeroslope=False
                slopeoffby = abs(slope+errorsig*slope_error)
                slopeoffby_rounded = round_to_1(slopeoffby)
        else:
            if (slope-errorsig*slope_error)<=0.0:
                zeroslope=True
            else:
                zeroslope=False
                slopeoffby = (slope-errorsig*slope_error)
                slopeoffby_rounded = round_to_1(slopeoffby)

    except RuntimeError:
        curvefit=False
        zeroslope=False
        #zeropoint,zeropoint_error  = float("nan"),float("nan")
        #slope,slope_error          = float("nan"),float("nan")
    
    if sc==False:
        newfname = file.replace(".fts","_magsoln_zeropoints.pdf")
    else:
        newfname = file.replace(".fts","_magsoln_zeropoints_sc.pdf")
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(mag_apass,magdiff,xerr=magerr_apass,yerr=magdiff_error,fmt="o",alpha=0.5,label="data")
    if plotmedian==True:
        plt.axhline(magdiff_median,color="black",linestyle="--",linewidth=2,label="median")
    plt.text(0.05,0.85,"zeropoint: "+str(magdiff_median_rounded)+" +/- "+str(magdiff_median_error_rounded),transform=ax.transAxes)
    if plotweightedavg==True:
        plt.axhline(magdiff_weightedavg,color="black",linestyle="-.",linewidth=2,label="weighted average")
    if curvefit==True:
        if yerrors==True:
            plt.plot(mag_apass_fit,deltamag_fit,color="purple",linewidth=3,label="curvefit with y errors")
        elif yerrors==False:
            plt.plot(mag_apass_fit,deltamag_fit,color="purple",linewidth=3,label="curvefit without y errors")
        plt.text(0.05,0.95,"slope: "+str(slope_rounded)+" +/- "+str(slope_error_rounded),transform=ax.transAxes)
        if zeroslope==False:
            plt.text(0.05,0.9,"slope off from zero by "+str(slopeoffby_rounded),transform=ax.transAxes)
    plt.xlabel("APASS Magnitude")
    plt.ylabel("(APASS Magnitude - SExtractor Magnitude)")
    plt.legend(prop={"size":10},loc="lower left",numpoints=1)
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)
    
    print "slope: "+str(slope_rounded)+" +/- "+str(slope_error_rounded)
    if zeroslope==True:
        print "slope is zero to within uncertainties"
    else:
        print "slope is not zero to wihin uncertainties; off by "+str(slopeoffby_rounded)
    #print "zeropoint: "+str(zeropoint_rounded)+" +/- "+str(zeropoint_error_rounded)
    print "zeropoint: "+str(magdiff_median_rounded)+" +/- "+str(magdiff_median_error_rounded)
    
    #return zeroslope,zeropoint,zeropoint_error
    return zeroslope,magdiff_median,magdiff_median_error

def fzeropoint_sc(ra,dec,x_world,y_world,mag_apass,magerr_apass,mag_sex,magerr_sex,clip_sigma):
    '''
    Measures the magnitude zero point of an image.
    '''
    zp            = mag_apass - mag_sex
    zp_sigmaclip  = sigma_clip(zp,sigma=clip_sigma,iters=5)
    zp_outliers   = zp_sigmaclip.mask
    zp_sc         = zp_sigmaclip[~zp_outliers].data
    
    # clean up using sigma-clipped mask
    ra_sc             = np.ma.masked_array(ra,mask=zp_outliers)[~zp_outliers].data
    dec_sc            = np.ma.masked_array(dec,mask=zp_outliers)[~zp_outliers].data
    x_world_sc        = np.ma.masked_array(x_world,mask=zp_outliers)[~zp_outliers].data
    y_world_sc        = np.ma.masked_array(y_world,mask=zp_outliers)[~zp_outliers].data
    mag_apass_sc      = np.ma.masked_array(mag_apass,mask=zp_outliers)[~zp_outliers].data
    magerr_apass_sc   = np.ma.masked_array(magerr_apass,mask=zp_outliers)[~zp_outliers].data
    mag_sex_sc        = np.ma.masked_array(mag_sex,mask=zp_outliers)[~zp_outliers].data
    magerr_sex_sc     = np.ma.masked_array(magerr_sex,mask=zp_outliers)[~zp_outliers].data

    return ra_sc,dec_sc,x_world_sc,y_world_sc,mag_apass_sc,magerr_apass_sc,mag_sex_sc,magerr_sex_sc,zp_sc

def plotmagsoln(file,mag_apass,magerr_apass,mag_sex,magerr_sex,sexplotfolderdir,sc=False):
    '''
    Plots the magnitude solution by comparing APASS to SExtractor:
    (SExtractor mag - APASS mag) vs APASS mag.
    '''
    
    magdiff,magdiff_error = fmagdiff(mag_sex,magerr_sex,mag_apass,magerr_apass)
    
    if sc==True:
        newfname = file.replace(".fts","_magsoln_sc.pdf")
    else:
        newfname = file.replace(".fts","_magsoln.pdf")
    
    # (APASS mag - SExtractor mag) vs. APASS mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(mag_apass,magdiff,xerr=magerr_apass,yerr=magdiff_error,fmt="o",alpha=0.5)
    plt.xlabel("APASS Magnitude")
    plt.ylabel("(APASS Magnitude - SExtractor Magnitude)")
    plt.savefig(sexplotfolderdir+newfname,bbox_inches="tight")
    plt.close(fig)

def sexmagsoln(file1,folder1,file2,folder2,plotdir,zeropoints,VERBOSE,linear=False):
    '''
    Plots the equivalent of (APASS mag - SExtractor mag) vs. APASS mag for two SExtractor images
    as well as the distribution of magnitude differences between the two catalogues.
    
    file1    = W3_1-S001-R001-C001-i_dupe-1-new-dsff.fts
    folder1  = 20111221
    file2    = W3_1-S001-R001-C001-i_dupe-2-new-dsff.fts
    folder2  = 20111221
    '''
    
    folderdir1 = findfolderdir(folder1,VERBOSE=VERBOSE)
    folderdir2 = findfolderdir(folder2,VERBOSE=VERBOSE)
    
    zp_zperr1 = zeropoints[folderdir1+file1][0]
    zp_zperr2 = zeropoints[folderdir2+file2][0]
    
    zp1,zperr1 = zp_zperr1[0],zp_zperr1[1]
    zp2,zperr2 = zp_zperr2[0],zp_zperr2[1]
    
    sexcat1 = file1.replace(".fts",".cat")
    sexcat2 = file2.replace(".fts",".cat")
    
    # SExtractor catalog data; file1
    sexdata1           = ascii.read(sexcatdir+folder1+"/"+sexcat1)
    x_world1           = sexdata1["X_WORLD"]               # barycenter position along world x-axis           [deg]
    y_world1           = sexdata1["Y_WORLD"]               # barycenter position along world y-axis           [deg]
    mag_iso1           = sexdata1["MAG_ISO"]               # isophotal manitude                               [mag]
    magerr_iso1        = sexdata1["MAGERR_ISO"]            # RMS uncertainty for ISO magnitude                [mag]
    flux_iso1          = sexdata1["FLUX_ISO"]              # flux density of ISO magnitude                    [ADU]
    fluxerr_iso1       = sexdata1["FLUX_ISO"]              # RMS uncertainty of ISO flux density              [ADU]
    mag_isocor1        = sexdata1["MAG_ISOCOR"]            # corrected isophotal magnitude                    [mag]
    magerr_isocor1     = sexdata1["MAGERR_ISOCOR"]         # RMS uncertainty for ISOCOR magnitude             [mag]
    flux_isocor1       = sexdata1["FLUX_ISOCOR"]           # flux density of ISOCOR magnitude                 [ADU]
    fluxerr_isocor1    = sexdata1["FLUX_ISOCOR"]           # RMS uncertainty of ISOCOR flux density           [ADU]
    mag_auto1          = sexdata1["MAG_AUTO"]              # kron-like elliptical aperture magnitude          [mag]
    magerr_auto1       = sexdata1["MAGERR_AUTO"]           # RMS uncertainty for AUTO magnitude               [mag]
    flux_auto1         = sexdata1["FLUX_AUTO"]             # flux density of AUTO magnitude                   [ADU]
    fluxerr_auto1      = sexdata1["FLUX_AUTO"]             # RMS uncertainty of AUTO flux density             [ADU]
    flags1             = sexdata1["FLAGS"]                 # extraction flags
    class_star1        = sexdata1["CLASS_STAR"]            # S/G classification
    
    # SExtractor catalog data; file2
    sexdata2           = ascii.read(sexcatdir+folder2+"/"+sexcat2)
    x_world2           = sexdata2["X_WORLD"]               # barycenter position along world x-axis           [deg]
    y_world2           = sexdata2["Y_WORLD"]               # barycenter position along world y-axis           [deg]
    mag_iso2           = sexdata2["MAG_ISO"]               # isophotal manitude                               [mag]
    magerr_iso2        = sexdata2["MAGERR_ISO"]            # RMS uncertainty for ISO magnitude                [mag]
    flux_iso2          = sexdata2["FLUX_ISO"]              # flux density of ISO magnitude                    [ADU]
    fluxerr_iso2       = sexdata2["FLUX_ISO"]              # RMS uncertainty of ISO flux density              [ADU]
    mag_isocor2        = sexdata2["MAG_ISOCOR"]            # corrected isophotal magnitude                    [mag]
    magerr_isocor2     = sexdata2["MAGERR_ISOCOR"]         # RMS uncertainty for ISOCOR magnitude             [mag]
    flux_isocor2       = sexdata2["FLUX_ISOCOR"]           # flux density of ISOCOR magnitude                 [ADU]
    fluxerr_isocor2    = sexdata2["FLUX_ISOCOR"]           # RMS uncertainty of ISOCOR flux density           [ADU]
    mag_auto2          = sexdata2["MAG_AUTO"]              # kron-like elliptical aperture magnitude          [mag]
    magerr_auto2       = sexdata2["MAGERR_AUTO"]           # RMS uncertainty for AUTO magnitude               [mag]
    flux_auto2         = sexdata2["FLUX_AUTO"]             # flux density of AUTO magnitude                   [ADU]
    fluxerr_auto2      = sexdata2["FLUX_AUTO"]             # RMS uncertainty of AUTO flux density             [ADU]
    flags2             = sexdata2["FLAGS"]                 # extraction flags
    class_star2        = sexdata2["CLASS_STAR"]            # S/G classification
    
    # obtain Source Extractor photometry type information
    if phottype=="ISO":
        # catalogue 1
        flux1        = flux_iso1
        fluxerr1     = fluxerr_iso1
        mag_sex1     = mag_iso1
        magerr_sex1  = magerr_iso1
        # catalogue 2
        flux2        = flux_iso2
        fluxerr2     = fluxerr_iso2
        mag_sex2     = mag_iso2
        magerr_sex2  = magerr_iso2
    elif phottype=="ISOCOR":
        # catalogue 1
        flux1        = flux_isocor1
        fluxerr1     = fluxerr_isocor1
        mag_sex1     = mag_isocor1
        magerr_sex1  = magerr_isocor1
        # catalogue 2
        flux2        = flux_isocor2
        fluxerr2     = fluxerr_isocor2
        mag_sex2     = mag_isocor2
        magerr_sex2  = magerr_isocor2
    elif phottype=="AUTO":
        # catalogue 1
        flux1        = flux_auto1
        fluxerr1     = fluxerr_auto1
        mag_sex1     = mag_auto1
        magerr_sex1  = magerr_auto1
        # catalogue 2
        flux2        = flux_auto2
        fluxerr2     = fluxerr_auto2
        mag_sex2     = mag_auto2
        magerr_sex2  = magerr_auto2

    # extract SExtractor R.A. and Dec.; file1
    sexra1        = x_world1
    sexdec1       = y_world1
    sexradec1     = [sexra1,sexdec1]
    sexradec1     = np.array(np.transpose(sexradec1))
    # extract SExtractor R.A. and Dec.; file2
    sexra2        = x_world2
    sexdec2       = y_world2
    sexradec2     = [sexra2,sexdec2]
    sexradec2     = np.array(np.transpose(sexradec2))
    # create KDTree and use it to match SExtractor catalogues
    sextree  = spatial.KDTree(sexradec2)         # create KD-Tree of catalogue 2
    matches  = sextree.query(sexradec1)          # find matching sources between catalogues 1 and 2
    dist     = np.array(matches[0])              # distance between nearest neighbours in degrees
    indices  = np.array(matches[1])              # indices that match catalogue 2 to catalogue 1
    
    # find matching sources
    mag_sex2_matches     = mag_sex2[indices]     # matches with 'mag_sex1'
    magerr_sex2_matches  = magerr_sex2[indices]  # matches with 'magerr_sex1'
    flux2_matches        = flux2[indices]        # matches with 'flux1'
    fluxerr2_matches     = fluxerr2[indices]     # matches with 'fluxerr1'
    flags2_matches       = flags2[indices]       # matches with 'flags1'
    
    # clean up photometry catalog
    conditions     = np.array( (dist<posmatch)     & 
                               (flux1>0.0)         & (fluxerr1>0.0)         & (flags1<=flagmax) &
                               (flux2_matches>0.0) & (fluxerr2_matches>0.0) & (flags2_matches<=flagmax) )
    i_clean        = np.array(np.where(conditions)[0])
    dist_clean     = dist[i_clean]
    indices_clean  = indices[i_clean]
    
    # clean catalogue 1
    mag_sex1_clean     = mag_sex1[i_clean]
    magerr_sex1_clean  = magerr_sex1[i_clean]
    flux1_clean        = flux1[i_clean]
    fluxerr1_clean     = fluxerr1[i_clean]
    flags1_clean       = flags1[i_clean]    
    # clean catalogue 2
    mag_sex2_clean     = mag_sex2_matches[i_clean]
    magerr_sex2_clean  = magerr_sex2_matches[i_clean]
    flux2_clean        = flux2_matches[i_clean]
    fluxerr2_clean     = fluxerr2_matches[i_clean]
    flags2_clean       = flags2_matches[i_clean]
    
    if VERBOSE:
        print len(sexradec1),     "stars extracted from SExtractor catalogue 1"
        print len(sexradec2),     "stars found using SExtractor catalogue 2"
        print len(indices_clean), "matching SExtractor stars"
    
    # correct for magnitude zero point
    mag_sex1_clean         = mag_sex1_clean + zp1
    mag_sex2_clean         = mag_sex2_clean + zp2
    magerr_sex1_clean      = np.sqrt( (magerr_sex1_clean)**2. + (zperr1)**2. )
    magerr_sex2_clean      = np.sqrt( (magerr_sex2_clean)**2. + (zperr2)**2. )
    
    magdiff,magdiff_error  = fmagdiff(mag_sex2_clean,magerr_sex2_clean,mag_sex1_clean,magerr_sex1_clean)
    sex_syserr             = (np.std(magdiff)/np.sqrt(len(magdiff)))/np.sqrt(2.)
    sex_syserr_rounded     = round_to_1(sex_syserr)
    
    if linear==True:
        fname = "sexmagdiff_linzpvsairmass.pdf"
    else:
        fname = "sexmagdiff_nonlinzpvsairmass.pdf"
    
    # SExtractor magnitude difference vs. calibrated SExtractor magnitude
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(mag_sex1_clean,magdiff,xerr=magerr_sex1_clean,yerr=magdiff_error,fmt="o",alpha=0.5)
    plt.text(0.05,0.95,"systematic error: "+str(sex_syserr_rounded)+" mag",transform=ax.transAxes)
    plt.text(0.05,0.9,"mean SExtractor magdiff error: "+str(np.mean(magdiff_error))+" mag",transform=ax.transAxes)
    plt.xlabel("Calibrated SExtractor Magnitude [mag]")
    plt.ylabel("Calibrated SExtractor Magnitude Difference")
    plt.savefig(plotdir+fname,bbox_inches="tight")
    plt.close(fig)
    
    # histogram of SExtractor magnitude difference
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    pylab.hist(magdiff,histtype="step",fill=True,edgecolor="black")
    plt.text(0.05,0.95,"systematic error: "+str(sex_syserr_rounded)+" mag",transform=ax.transAxes)
    plt.xlabel("SExtractor Magnitude Difference")
    plt.ylabel("N")
    plt.savefig(plotdir+"sexmagdiff_hist.pdf",bbox_inches="tight")
    plt.close(fig)
    
def zpvsairmass(zeropoints,plotdir):
    '''
    Plots the magnitude zero point versus airmass.
    
    fig1 : zero point versus airmass
    fig2 : zero point uncertainty versus airmass
    '''
    
    zp_zperr_airmass_dict = {}
    passbands             = []
    
    for filedir in zeropoints.keys():
        zp_zperr  = zeropoints[filedir][0]
        zp,zperr  = zp_zperr[0],zp_zperr[1]
        header    = getheader(filedir)
        airmass   = header["AIRMASS"]
        passband  = header["FILTER"]
        if passband not in passbands:
            passbands.append(passband)
        appenditemtodict(zp,"zp_"+passband,zp_zperr_airmass_dict)
        appenditemtodict(zperr,"zperr_"+passband,zp_zperr_airmass_dict)
        appenditemtodict(airmass,"airmass_"+passband,zp_zperr_airmass_dict)
        
        if passband=="i":
            if airmass>1.22:
                print "linear airmass relation:", filedir
            elif airmass<1.15:
                print "non-linear airmass relation:", filedir
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212,sharex=ax1)
    
    for passband in passbands:
        zp_all       = zp_zperr_airmass_dict["zp_"+passband]
        zperr_all    = zp_zperr_airmass_dict["zperr_"+passband]
        airmass_all  = zp_zperr_airmass_dict["airmass_"+passband]
        ax1.errorbar(airmass_all,zp_all,yerr=zperr_all,fmt="o",label=passband+"-band")
        ax2.errorbar(airmass_all,zperr_all,fmt="o")
    
    # hide x-ticks on top subplot since both subplots share the x-axis
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    ax2.set_xlabel("Airmass (fraction of zenith airmass)")
    ax1.set_ylabel("Zero Point [mag]")
    ax2.set_ylabel("Uncertainty in Zero Point [mag]")
    ax1.legend(loc="lower right",numpoints=1)
    plt.subplots_adjust(hspace=0.075)
    plt.savefig(plotdir+"zpvsairmass.pdf",bbox_inches="tight")
    plt.close(fig)
    
def pointingoffsethist(file,apassra,apassdec,sexra,sexdec,indices,dist,sexplotfolderdir,radius=posmatch,deltapos=1.0,poslower=0.0,posupper=10.0):
    '''
    Plots the pointing offset between APASS and 
    Source Extractor catalog matches from KDTree.
    '''
    
    apassra_matches   = apassra[indices]
    apassdec_matches  = apassdec[indices]
    
    deltara           = abs(np.array(apassra_matches-sexra))
    deltadec          = abs(np.array(apassdec_matches-sexdec))
    deltapos_deg      = np.sqrt((deltara)**2. + (deltadec)**2.)
    deltapos_arcsec   = deltapos_deg*3600.
    dist_arcsec       = dist*3600.
    
    newfname = file.replace(".fts","_posresidualhist.pdf")
    
    N = len(np.arange(poslower,posupper,deltapos))
    params = dict(bins=N,
                  range=(poslower,posupper))
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    pylab.hist(deltapos_arcsec,**params)
    plt.axvline(radius,linestyle="dashed",color="black",label="position match")
    plt.xlabel("Residual Pointing Offset (arcseconds)")
    plt.ylabel("N")
    plt.legend(loc="upper right")
    plt.savefig(sexplotfolderdir+newfname,bbox_inches="tight")
    plt.close(fig)
