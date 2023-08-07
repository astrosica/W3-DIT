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

parentdir       = config_data["parentdir"]
plotdir         = parentdir+"plots/"
sexplotdir      = plotdir+"sextractor/"
phottype        = overlapconfig_data["phottype"]                     # SExtractor photometry type; one of ISO, ISOCOR, or AUTO
detect_thresh   = float(overlapconfig_data["detect_thresh"])         # SNR threshold
posmatch        = float(overlapconfig_data["posmatch"])/3600.        # convert to degrees
posmatch_lower  = float(overlapconfig_data["posmatch_lower"])/3600.  # convert to degrees
posmatch_upper  = float(overlapconfig_data["posmatch_upper"])/3600.  # convert to degrees
magmin          = float(overlapconfig_data["magmin"])                # minimum APASS magnitude
magmax          = float(overlapconfig_data["magmax"])                # maximum APASS magnitude
flagmax         = float(overlapconfig_data["flagmax"])               # maximum Source Extractor internal flag
clip_sigma      = float(overlapconfig_data["clip_sigma"])            # # of standard deviations for upper and lower clipping limits

#########################################################################################################
#                                                                                                       #
#                                         IMPORT PACKAGES                                               #
#                                                                                                       #
#########################################################################################################


import os
import numpy as np
from scipy import spatial
from astropy.io import ascii
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
#                                       SOURCE EXTRACTOR                                                #
#                                      CONFIGURATION FILES                                              #
#                                                                                                       #
#########################################################################################################


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
DETECT_THRESH    {detect_thresh}              # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
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

default_conv = '''CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
'''

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


def fgauss(x, amp, mean, sigma):
    '''
    Gaussian distribution function.
    x      : input variable
    amp    : maximum amplitude
    mean   : distribution mean
    sigma  : standard deviation
    '''
    y = amp * np.exp(-((x - mean)**2. / (2. * sigma**2.)))
    return y


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
	f	 : file name to run source extractor on
	fpath	 : directory path to f
	objdir   : directory to save object map in
	catdir   : directory to save catalogue in
	bgdir	 : directory to save background map in
	Returns location of output catalog
	'''
	
        # Split to make file name for catalogue, object map and background map filenames
	fname = f.split(".fts")[0]
        if sexbool==True:
            # Construct source extractor call
            objsexcall = "sex -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME "+objdir+folder+"/"+fname+"_objects.fts "+fpath+f
            baksexcall = "sex -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME default.param -CATALOG_NAME "+catdir+folder+"/"+fname+".cat"+" -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME "+bgdir+folder+"/"+fname+"_background.fts "+fpath+f
            os.system(objsexcall)
            os.system(baksexcall)
        
        return catdir+folder+"/"+fname+".cat"

def apassoverlap(file,folder,folderdir,apassdata,sexcat,overlapannfolderdir,VERBOSE,GENERATE):
    '''
    Builds KDTree of APASS catalog to be compared with Source Extractor catalogue.
    '''
    
    # load SExtractor catalog data
    sexdata         = ascii.read(sexcat)
    x_world         = sexdata["X_WORLD"]
    errx2_world     = sexdata["ERRX2_WORLD"]
    y_world         = sexdata["Y_WORLD"]
    erry2_world     = sexdata["ERRY2_WORLD"]
    
    # extract SExtractor R.A. and Dec.
    sexra           = x_world
    sexdec          = y_world
    sexradec        = [sexra,sexdec]
    sexradec        = np.array(np.transpose(sexradec))
    # extract APASS R.A. and Dec.
    apassra          = apassdata[:,0]
    apassdec         = apassdata[:,2]
    apassradec       = [apassra,apassdec]
    apassradec       = np.array(np.transpose(apassradec))
    
    # create APASS KDTree and use it to query SExtractor sources
    apasstree   = spatial.KDTree(apassradec)  # create KD-Tree
    matches     = apasstree.query(sexradec)   # find which APASS sources are closest to each DIT source
    dist        = np.array(matches[0])        # distance between nearest neighbours in degrees
    indices     = np.array(matches[1])        # APASS indices that match with SExtractor sources
    
    if VERBOSE:
        print len(apassra), "stars extracted from APASS catalog"
        print len(dist), "matching stars extracted from SExtractor"
    
    if GENERATE:
        sexplotfolderdir = sexplotdir+folder+"/"
        # histogram of pointing offsets between APASS and SExtractor catalogs
        pointingoffsethist(file,apassra,apassdec,sexra,sexdec,indices,dist,sexplotfolderdir)

    # clean up APASS overlap
    fphotmatches(file,folder,folderdir,apassdata,indices,dist,sexcat,sexplotfolderdir,overlapannfolderdir,VERBOSE,GENERATE)

def fphotmatches(
    file,folder,folderdir,apassdata,indices,dist,sexcat,sexplotfolderdir,overlapannfolderdir,VERBOSE,GENERATE,magmin=magmin,magmax=magmax,flagmax=flagmax):
    '''
    Cleans up APASS and Source Extractor photometry data.
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
    #snr_iso           = abs(flux_iso/fluxerr_iso)        # ISO flux density SNR
    #snr_isocor        = abs(flux_isocor/fluxerr_isocor)  # ISOCOR flux density SNR
    #snr_auto          = abs(flux_auto/fluxerr_auto)      # AUTO flux density SNR
    
    header            = getheader(folderdir+file)
    passband          = header["FILTER"]
     
    # APASS catalog data
    ra                = apassdata[:,0]     # right ascension in degrees
    ra_err            = apassdata[:,1]     # error in right ascension in arcseconds
    dec               = apassdata[:,2]     # declination in degrees
    dec_err           = apassdata[:,3]     # error in declination in arcseconds
    Nobs              = apassdata[:,4]     # number of observations
    V_mag             = apassdata[:,5]     # Johnson V magnitude
    V_magerr          = apassdata[:,6]     # error in Johnson V magnitude
    B_mag             = apassdata[:,7]     # Johnson B magnitude
    B_magerr          = apassdata[:,8]     # error in Johnson B magnitude
    g_mag             = apassdata[:,9]     # Sloan g' magnitude
    g_magerr          = apassdata[:,10]    # error in Sloan g' magnitude
    r_mag             = apassdata[:,11]    # Sloan r' magnitude
    r_magerr          = apassdata[:,12]    # error in Sloan r' magnitude
    i_mag             = apassdata[:,13]    # Sloan i' magnitude
    i_magerr          = apassdata[:,14]    # error in Sloan i' magnitude
    
    # APASS matches with SExtractor sources
    ra_matches        = ra[indices]        # right ascension in degrees
    ra_err_matches    = ra_err[indices]    # error in right ascension in arcseconds
    dec_matches       = dec[indices]       # declination in degrees
    dec_err_matches   = dec_err[indices]   # error in declination in arcseconds
    Nobs_matches      = Nobs[indices]      # number of observations
    V_mag_matches     = V_mag[indices]     # Johnson V magnitude
    V_magerr_matches  = V_magerr[indices]  # error in Johnson V magnitude
    B_mag_matches     = B_mag[indices]     # Johnson B magnitude
    B_magerr_matches  = B_magerr[indices]  # error in Johnson B magnitude
    g_mag_matches     = g_mag[indices]     # Sloan g' magnitude
    g_magerr_matches  = g_magerr[indices]  # error in Sloan g' magnitude
    r_mag_matches     = r_mag[indices]     # Sloan r' magnitude
    r_magerr_matches  = r_magerr[indices]  # error in Sloan r' magnitude
    i_mag_matches     = i_mag[indices]     # Sloan i' magnitude
    i_magerr_matches  = i_magerr[indices]  # error in Sloan i' magnitude
    
    # obtain passband information
    if passband=="i":
        mags     = i_mag_matches
        magerrs  = i_magerr_matches
    elif passband=="r":
        mags     = r_mag_matches
        magerrs  = r_magerr_matches
    elif passband=="z":
        mags     = g_mag_matches
        magerrs  = g_magerr_matches

    # obtain Source Extractor photometry type information
    if phottype=="ISO":
        flux = flux_iso
        fluxerr = fluxerr_iso
        
    elif phottype=="ISOCOR":
        flux = flux_isocor
        fluxerr = fluxerr_isocor
        
    elif phottype=="AUTO":
        flux = flux_auto
        fluxerr = fluxerr_auto

    
    # clean up photometry catalog
    conditions = np.array(
        (dist<posmatch) & (flux>0.0) & (fluxerr>0.0) & (mags<=magmax) & (mags>=magmin) & (flags<=flagmax)
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
    mags_clean            = mags[i]
    magerrs_clean         = magerrs[i]
    
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
    #snr_iso_clean         = snr_iso[i]
    #snr_isocor_clean      = snr_isocor[i]
    #snr_auto_clean        = snr_auto[i]

    # perform sigma-clipping on magnitude zero points
    ra_sc,dec_sc,mags_sc,magerrs_sc,mag_isocor_sc,magerr_isocor_sc,zp_sc=fzeropoint(ra_clean,dec_clean,mags_clean,magerrs_clean,mag_isocor_clean,magerr_isocor_clean,clip_sigma)
    
    # perform linear fitting to calculate magnitude zero point and slope
    '''
    fitzeropoint(mags_sc,magerrs_sc,mag_isocor_sc,magerr_isocor_sc)
    '''
    
    if GENERATE:
        # plot histograms of (APASS mag - SExtractor mag)
        plotmagsolnhist(file,mags_clean,mag_isocor_clean,dist_clean,sexplotfolderdir)
        
        # plot histograms of APASS and SExtractor mag uncertainties
        plotmagerrorhist(file,magerrs_clean,magerr_iso_clean,magerr_isocor_clean,magerr_auto_clean,sexplotfolderdir)

        # plot (APASS mag- SExtractor mag) vs. APASS mag
        plotmagsolution(file,mags_clean,magerrs_clean,mag_isocor_clean,magerr_isocor_clean,sexplotfolderdir)
        #plotmagsolution(file,mags_sc,magerrs_sc,mag_isocor_sc,magerr_isocor_sc,sexplotfolderdir,sc=True)
        
        # plot ISO, ISOCOR, and AUTO magnitudes
        comparemags(file,mag_iso_clean,magerr_iso_clean,mag_isocor_clean,magerr_isocor_clean,mag_auto_clean,magerr_auto_clean,mags_clean,magerrs_clean,sexplotfolderdir)
        
        # plot SExtractor and APASS positional errors versus SExtractor's ISO, ISOCOR, and AUTO magnitudes
        poserrvsmags(file,ra_err_clean,dec_err_clean,xerr_world_clean,yerr_world_clean,mag_iso_clean,magerr_iso_clean,mag_isocor_clean,magerr_isocor_clean,mag_auto_clean,magerr_auto_clean,mags_clean,magerrs_clean,sexplotfolderdir)
        
        # create annotation files for APASS matches
        writeannfiles(file,ra_clean,dec_clean,"_clean",overlapannfolderdir,VERBOSE)   # cleaned
        #writeannfiles(file,ra_sc,dec_sc,"_clean_sc",overlapannfolderdir,VERBOSE)      # cleaned and sigma-clipped

def plotmagerrorhist(file,magerrs,magerr_iso,magerr_isocor,magerr_auto,plotdir):
    '''
    Plot histograms of APASS and SExtractor magnitude uncertainties.
    '''
    
    deltaerr = 0.05
    errlower = np.min([np.min(magerrs),np.min(magerr_iso),np.min(magerr_isocor),np.min(magerr_auto)])
    errupper = np.max([np.max(magerrs),np.max(magerr_iso),np.max(magerr_isocor),np.max(magerr_auto)])
    
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
    pylab.hist((magerrs,magerr_iso,magerr_isocor,magerr_auto),**commonparams)
    plt.xlabel("Magnitude Uncertainty")
    plt.xlim(errlower,errupper)
    plt.ylabel("N")
    plt.legend(loc="upper right",handles=[legend_apass,legend_iso,legend_isocor,legend_auto],prop={"size":14})
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close()
    plt.clf()

def plotmagsolnhist(file,mags,mag_isocor,dist,plotdir,deltamag=0.05):
    '''
    Plots histograms of (APASS mag - SExtractor mag) for matches
    within both 1 and 5 arcseconds.
    '''
    
    condition_1arcsec    = np.array(dist<=1.0/3600.)
    condition_5arcsec    = np.array( (dist>1.0/3600.) & (dist<=5.0/3600.) )
    i_1arcsec            = np.array(np.where(condition_1arcsec)[0])
    i_5arcsec            = np.array(np.where(condition_5arcsec)[0])
    mags_1arcsec         = mags[i_1arcsec]
    mags_5arcsec         = mags[i_5arcsec]
    mag_isocor_1arcsec   = mag_isocor[i_1arcsec]
    mag_isocor_5arcsec   = mag_isocor[i_5arcsec]
    magdiff_1arcsec      = mags_1arcsec - mag_isocor_1arcsec
    magdiff_5arcsec      = mags_5arcsec - mag_isocor_5arcsec
    
    #deltamag = 0.05
    maglower = np.min([np.min(magdiff_1arcsec),np.min(magdiff_1arcsec)])
    magupper = np.max([np.max(magdiff_5arcsec),np.max(magdiff_5arcsec)])
    
    histmin = np.min([np.min(magdiff_1arcsec),np.min(magdiff_5arcsec)])
    histmax = np.max([np.max(magdiff_1arcsec),np.max(magdiff_5arcsec)])

    N = len(np.arange(histmin,histmax,deltamag))
    commonparams = dict(bins=N,range=(histmin,histmax))
    histbins = np.linspace(histmin,histmax,N,endpoint=True)
    
    # fit Gaussian
    
    yhist_1arcsec,xhist_1arcsec = np.histogram(magdiff_1arcsec,bins=histbins)
    yhist_5arcsec,xhist_5arcsec = np.histogram(magdiff_5arcsec,bins=histbins)
    
    xhist_halfwidth_1arcsec = 0.5*abs(xhist_1arcsec[1]-xhist_1arcsec[0])
    xhist_halfwidth_5arcsec = 0.5*abs(xhist_5arcsec[1]-xhist_5arcsec[0])
    
    xhist_gauss_1arcsec = xhist_1arcsec[:-1]+xhist_halfwidth_1arcsec
    xhist_gauss_5arcsec = xhist_5arcsec[:-1]+xhist_halfwidth_5arcsec
    
    mu_1arcsec, sigma_1arcsec = norm.fit(magdiff_1arcsec)
    mu_5arcsec, sigma_5arcsec = norm.fit(magdiff_5arcsec)
    
    try:
        fit=True
        popt_1arcsec,pcov_1arcsec = curve_fit(
            fgauss,xhist_gauss_1arcsec,yhist_1arcsec,np.array([np.max(yhist_1arcsec),mu_1arcsec,sigma_1arcsec]),maxfev=1000)
        popt_5arcsec,pcov_5arcsec = curve_fit(
            fgauss,xhist_gauss_5arcsec,yhist_5arcsec,np.array([np.max(yhist_5arcsec),mu_5arcsec,sigma_5arcsec]),maxfev=1000)
    except RuntimeError:
        fit=False
    
    newfname = file.replace(".fts","_magsoln_hist.pdf")
    
    legend_1arcsec = mpatches.Patch(color="blue",label="1 arcsec")
    legend_5arcsec = mpatches.Patch(color="green",label="5 arcsec")

    # histogram of (APASS mag - SExtractor mag)
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    # histograms
    pylab.hist((magdiff_1arcsec,magdiff_5arcsec),**commonparams)
    # fits
    if fit==True:
        plt.plot(np.linspace(histmin,histmax,100),fgauss(np.linspace(histmin,histmax,100),*popt_1arcsec),"b--",linewidth=2)
        plt.plot(np.linspace(histmin,histmax,100),fgauss(np.linspace(histmin,histmax,100),*popt_5arcsec),"g--",linewidth=2)
        plt.axvline(popt_1arcsec[1],color="blue",linestyle="dashed",linewidth=2)
        plt.axvline(popt_5arcsec[1],color="green",linestyle="dashed",linewidth=2)
    elif fit==False:
        pass
    plt.xlabel("(APASS Magnitude - SExtractor Magnitude)")
    plt.xlim(histmin,histmax)
    plt.ylabel("N")
    plt.legend(loc="upper left",handles=[legend_1arcsec,legend_5arcsec],prop={"size":14})
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close()
    plt.clf()
    
def poserrvsmags(file,ra_err,dec_err,xerr_world,yerr_world,mag_iso,magerr_iso,mag_isocor,magerr_isocor,mag_auto,magerr_auto,mags,magerrs,plotdir):
    '''
    Plots Source Extractor's positional error versus APASS magnitudes.
    '''
    
    poserr_APASS = np.sqrt((ra_err)**2. + (dec_err)**2.)
    poserr_SEX   = np.sqrt((xerr_world)**2. + (yerr_world)**2.)*3600.
    
    newfname_SEX   = file.replace(".fts","_poserrvsmags_sextractor.pdf")
    newfname_APASS = file.replace(".fts","_poserrvsmags_apass.pdf")
    
    fig = plt.figure(figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(mag_iso,poserr_SEX,xerr=magerr_iso,fmt="o",color="purple",label="ISO")
    plt.errorbar(mag_isocor,poserr_SEX,xerr=magerr_isocor,fmt="o",color="orange",label="ISOCOR")
    plt.errorbar(mag_auto,poserr_SEX,xerr=magerr_auto,fmt="o",color="turquoise",label="AUTO")
    # plot log
    ax.set_yscale("log")
    plt.xlabel("Source Extractor Magnitude")
    plt.ylabel("Sextractor RMS Positional Uncertainty")
    plt.legend(loc="upper left",numpoints=1)
    plt.savefig(plotdir+newfname_SEX,bbox_inches="tight")
    plt.close()
    plt.clf()
    
    fig = plt.figure(figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(mags,poserr_APASS,xerr=magerrs,fmt="o",color="purple")
    # plot log
    ax.set_yscale("log")
    plt.xlabel("APASS Magnitude")
    plt.ylabel("APASS RMS Positional Uncertainty")
    plt.savefig(plotdir+newfname_APASS,bbox_inches="tight")
    plt.close()
    plt.clf()

def comparemags(file,mag_iso,magerr_iso,mag_isocor,magerr_isocor,mag_auto,magerr_auto,mags,magerrs,plotdir):
    '''
    Plots Source Extractor's ISO, ISOCOR, and AUTO magnitudes versus APASS magnitudes.
    '''
    
    # APASS mag - SExtractor mag
    magdiff_iso           = mags - mag_iso
    magdiff_isocor        = mags - mag_isocor
    magdiff_auto          = mags - mag_auto
    # uncertainties in APASS mag - SExtractor mag
    magdiff_iso_error     = np.sqrt((magerrs)**2. + (magerr_iso)**2.)
    magdiff_isocor_error  = np.sqrt((magerrs)**2. + (magerr_isocor)**2.)
    magdiff_auto_error    = np.sqrt((magerrs)**2. + (magerr_auto)**2.)
    
    newfname = file.replace(".fts","_comparemags.pdf")
    
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    # (APASS mag - SExtractor mag) vs. APASS mag
    plt.errorbar(mags,magdiff_iso,xerr=magerrs,yerr=magdiff_iso_error,fmt="o",color="purple",label="ISO")
    plt.errorbar(mags,magdiff_isocor,xerr=magerrs,yerr=magdiff_isocor_error,fmt="o",color="orange",label="ISOCOR")
    plt.errorbar(mags,magdiff_auto,xerr=magerrs,yerr=magdiff_auto_error,fmt="o",color="turquoise",label="AUTO")
    '''
    # SExtractor mag vs. APASS mag
    plt.errorbar(mags,mag_iso,xerr=magerrs,yerr=magerr_iso,fmt="o",color="purple",label="ISO")
    plt.errorbar(mags,mag_isocor,xerr=magerrs,yerr=magerr_isocor,fmt="o",color="orange",label="ISOCOR")
    plt.errorbar(mags,mag_auto,xerr=magerrs,yerr=magerr_auto,fmt="o",color="turquoise",label="AUTO")
    '''
    plt.xlabel("APASS Magnitude")
    plt.ylabel("Source Extractor Magnitude")
    plt.legend(loc="upper left",numpoints=1)
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close()
    plt.clf()

def apasshist(apassdata,maglower,magupper,deltamag,plotdir):
    '''
    Create a histogram of APASS magnitudes.
    
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
        mags          = apassdata[:,magcol]
        magerrs       = apassdata[:,magerrcol]
        magstring     = magstrings[i]
        magerr_fname  = magerr_fnames[i]
        colour        = colours[i]
        
        plt.plot(mags,magerrs,"ro",color=colour,alpha=0.5,label=magstring)
    
    plt.axhline(0,linestyle="dashed",color="black")
    plt.xlabel("Magnitude")
    plt.ylabel("Magnitude Uncertainty")
    plt.legend(loc="upper left",numpoints=1)
    plt.savefig(plotdir+"allmagerr.pdf",bbox_inches="tight")
    plt.close()
    plt.clf()
    
    # combined magnitude histogram plot using input limits
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    
    hist1 = apassdata[:,magcols[0]]
    hist2 = apassdata[:,magcols[1]]
    hist3 = apassdata[:,magcols[2]]
    hist4 = apassdata[:,magcols[3]]
    hist5 = apassdata[:,magcols[4]]
    
    N = len(np.arange(maglower,magupper,deltamag))
    commonparams = dict(bins=N,
                        range=(maglower,magupper))
    
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
    plt.clf()
    
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
    commonparams = dict(bins=N,
                        range=(maglower,magupper))
    
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
    plt.close()
    plt.clf()

def writeannfiles(file,ra_all,dec_all,fname_suffix,overlapannfolderdir,VERBOSE,colour="RED"):
    '''
    Writes APASS matches to annotation files in both
    kvis and ds9 formats.
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

def fitzeropoint(mags_sc,magerrs_sc,mag_isocor_sc,magerr_isocor_sc):
    '''
    Fits a straight line to sigma-clipped(SExtractor - APASS) vs. APASS magnitudes
    to measure magnitude zero point and slope in the data.
    '''
    
    def linearfunc(m,x,b):
        return m*x + b
    

def fzeropoint(ra_clean,dec_clean,mags_clean,magerrs_clean,mag_isocor_clean,magerr_isocor_clean,clip_sigma):
    '''
    Measures the magnitude zero point of an image.
    '''
    zp_clean      = mag_isocor_clean - mags_clean
    zp_sigmaclip  = sigma_clip(zp_clean,sigma=clip_sigma,iters=5)
    zp_outliers   = zp_sigmaclip.mask
    zp_sc         = zp_sigmaclip[~zp_outliers].data
    
    # clean up using sigma-clipped mask
    ra_sc             = np.ma.masked_array(ra_clean,mask=zp_outliers)[~zp_outliers].data
    dec_sc            = np.ma.masked_array(dec_clean,mask=zp_outliers)[~zp_outliers].data
    mags_sc           = np.ma.masked_array(mags_clean,mask=zp_outliers)[~zp_outliers].data
    magerrs_sc        = np.ma.masked_array(magerrs_clean,mask=zp_outliers)[~zp_outliers].data
    mag_isocor_sc     = np.ma.masked_array(mag_isocor_clean,mask=zp_outliers)[~zp_outliers].data
    magerr_isocor_sc  = np.ma.masked_array(magerr_isocor_clean,mask=zp_outliers)[~zp_outliers].data

    return ra_sc,dec_sc,mags_sc,magerrs_sc,mag_isocor_sc,magerr_isocor_sc,zp_sc

def plotmagsolution(file,mags,magerrs,mag_isocor,magerr_isocor,sexplotfolderdir,sc=False):
    '''
    Plots the magnitude solution by comparing APASS to SExtractor:
    (SExtractor mag - APASS mag) vs APASS mag.
    '''
    
    apassmag        = mags
    apassmag_error  = magerrs
    sexmag          = mag_isocor
    sexmag_error    = magerr_isocor
    
    magdiff         = apassmag-sexmag
    magdiff_error   = np.sqrt((sexmag_error)**2. + (apassmag_error)**2.)
    
    if sc==True:
        newfname = file.replace(".fts","_magsoln_sc.pdf")
    else:
        newfname = file.replace(".fts","_magsoln.pdf")
    
    Npoints = len(magdiff)
    
    # (APASS mag - SExtractor mag) vs. APASS mag
    fig = plt.figure(1,figsize=(11,8.5))
    ax = fig.add_subplot(111)
    plt.errorbar(apassmag,magdiff,xerr=apassmag_error,yerr=magdiff_error,fmt="o",alpha=0.5)
    plt.xlabel("APASS Magnitude")
    plt.ylabel("(APASS Magnitude - SExtractor Magnitude)")
    plt.title(str(Npoints))
    plt.savefig(sexplotfolderdir+newfname,bbox_inches="tight")
    plt.close()
    plt.clf()
    
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
    plt.close()
    plt.clf()
