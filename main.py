#!/usr/local/bin/python

'''
Requires that images are organized in the following way:
raw science images: <parent directory>/<date>/<images>/
bias images: <parent directory>/<date>/<calibration folder>/<images>/
dark images: <parent directory>/<date>/<calibration folder>/<images>/
flat images: <parent directory>/<date>/<flat folder>/<images>/

Requires that filenames are named in the following way:
raw science images: objectname_[]-[]-[].fts
bias images: Bias-[]-[]-[]-.fts
dark images: Dark-[]-[]-[].fts
flat images: AutoFlat-[]-[]-[].fts

The following is done with the raw science images:
1. dark subtracted and flat fielded
2. astronometry.net is used to obtain accurate world coordinate system
3. Source Extractor is used for photometry

Outputs the following directories
- :

Requires the following header keys in the science images' FITS files:
- OBJECT
- FILTER
- EXPTIME

Requires the following software:
    - sextractor
    - astrometry.net
Requires the following packages: os, re, yaml, docopt, numpy, astropy
Requires the following files:
    - calibration.py
    - callastrometry.py

Usage:
main [options]

Options:
  -h, --help           Show this screen.
                       [default: False]
  -l, --logoutput      Generate log outputs.
                       [default: True]
  -v, --verbose        Generate verbose output.
                       [default: True]
  -g, --generate       Generate outputs.
                       [default: True]
  -w, --wcs            Create WCS solution via Astronomy.net.
                       [default: False]
  -r, --reduce         Perform data reduction.
                       [default: False]
  -c, --calibrate      Perform zero-point calibration.
                       [default: False]
  -p, --panstarrs      Use Pan-STARRS for calibration.
                       [default: False]
  -a, --atlas          Use ATLAS for calibration.
                       [default: True]
  -s, --sextractor     Run Source Extractor for calibration step.
                       [default: False]
  -k, --checkcal       Check the zero-point calibration.
                       [default: False]
  -x, --stack          Stack calibrated images.
                       [default: False]
  -f, --finalcalcheck  Check the zero-point calibration on stacked images.
                       [default: False]
  -n, --finalsexrun    Final SExtractor run on stacked images.
                       [default: False]
  -m, --matchcats      Match source catalogues.
                       [default: False]
'''


#########################################################################################################
#                                                                                                       #
#                                         IMPORT PACKAGES                                               #
#                                                                                                       #
#########################################################################################################


import os
import re
import yaml
import docopt
import cPickle as pickle
from astropy.io import fits
from astropy.table import Table
import subprocess


#########################################################################################################
#                                                                                                       #
#                                    COMMAND LINE ARGUMENTS                                             #
#                                                                                                       #
#########################################################################################################

arguments     = docopt.docopt(__doc__)
# options without arguments
HELP          = arguments["--help"]
LOGOUTPUT     = arguments["--logoutput"]
VERBOSE       = arguments["--verbose"]
GENERATE      = arguments["--generate"]
CREATEWCS     = arguments["--wcs"]
REDUCE        = arguments["--reduce"]
CALIBRATE     = arguments["--calibrate"]
PANSTARRS     = arguments["--panstarrs"]
ATLAS         = arguments["--atlas"]
SEXTRACTOR    = arguments["--sextractor"]
CHECKCAL      = arguments["--checkcal"]
STACK         = arguments["--stack"]
FINALCALCHECK = arguments["--finalcalcheck"]
FINALSEXRUN   = arguments["--finalsexrun"]
MATCHCATS     = arguments["--matchcats"]

#########################################################################################################
#                                                                                                       #
#                                 IMPORT CONFIGURATION FILE                                             #
#                                                                                                       #
#########################################################################################################


with open("config.yml","r") as fconfig:
    config_data = yaml.safe_load(fconfig)

# input:
ncores_wcs               = int(config_data["ncores_wcs"])          # number of cores for WCS solution if multiprocessing is used
ncores_cal               = int(config_data["ncores_cal"])          # number of cores for data reduction if multiprocessing is used
ncores_apass             = int(config_data["ncores_apass"])        # number of cores for APASS solution if multiprocessing is used
ncores_panstarrs         = int(config_data["ncores_panstarrs"])    # number of cores for Pan-STARRS solution if multiprocessing is used
ncores_atlas             = int(config_data["ncores_atlas"])        # number of cores for ATLAS solution if multiprocessing is used
MULTIPROCESSING_CAL      = config_data["MULTIPROCESSING_CAL"]      # multiprocessing boolean for calibration (if True, excute multiprocessing)

parentdir                = config_data["parentdir"]                # parent directory where code is stored
datadir                  = config_data["datadir"]                  # directory where data is stored
calfolder                = config_data["calfolder"]                # folder name containing bias and dark images
flatfolder               = config_data["flatfolder"]               # folder name containing flat images
objects                  = config_data["objects"].split(" ")       # list of telescope pointings/object names
passbands                = config_data["passbands"].split(" ")     # list of passbands to be used
sciencedata_dsff_cutoff  = config_data["sciencedata_dsff_cutoff"]  # lower-limit cut-off value for dark-subtracted, flat-fielded science images
totalFOV                 = config_data["totalFOV"]                 # total field-of-view (FOV) covering all telescope pointings                  [deg]
centerRA                 = config_data["centerRA"]                 # Right Ascension (RA) of the center of all telescope pointings               [deg]
centerDEC                = config_data["centerDEC"]                # Declination (Dec) of the center of all telescope pointings                  [deg]
# output:
dsffdir                  = config_data["dsffdir"]                  # directory where reduced dark-subtracted, flat-fielded science images will be stored
mastercaldir             = config_data["mastercaldir"]             # directory where master 'calibration' images will be stored for data reduction
MMcaldir                 = mastercaldir+"MMcalimages/"             # directory where master-master 'calibration' images will be stored for data
calibrateddir            = config_data["calibrateddir"]            # directory where calibrated images will be stored
stackedcaldir            = config_data["stackedcaldir"]            # directory where final stacked calibrated images will be stored
# reduction
sciencemaskdir           = config_data["sciencemaskdir"]           # directory where science masks will be stored for data reduction
wcsdir                   = config_data["wcsdir"]                   # directory where Astrometry.net WCS solutions will be stored
apassdir                 = config_data["apassdir"]                 # directory where APASS catalogues will be stored
panstarrsdir             = config_data["panstarrsdir"]             # directory where PanSTARRS catalogue is stored
panstarrscat             = config_data["panstarrscat"]             # PanSTARRS catalogue file name
panstarrstrimcat         = config_data["panstarrstrimcat"]         # Trimmed PanSTARRS catalogue file name
panstarrsphot            = config_data["panstarrsphot"]            # Pan-STARRS photometry type: "Ap" (aperture), "PSF", or "Kron"
atlasdir                 = config_data["atlasdir"]                 # directory where ATLAS catalogue is stored
atlascat                 = config_data["atlascat"]                 # ATLAS catalogue file name
atlastrimcat             = config_data["atlastrimcat"]             # ATLAS catalogue file name
gaiadir                  = config_data["gaiadir"]                  # directory where Gaia catalogue is stored
gaiacat                  = config_data["gaiacat"]                  # Gaia catalogue file name
mergedcatsdir            = config_data["mergedcatsdir"]            # directory where merged catalogues are stored
ATLASDITcat              = config_data["ATLASDITcat"]              # ATLAS-DIT merged catalogue file name

# SExtractor directories
sexdir                   = config_data["sexdir"]                   # directory where Source Extractor files will be stored
sexobjdir                = config_data["sexobjdir"]                # directory where Source Extractor object files will be stored
sexcatdir                = config_data["sexcatdir"]                # directory where Source Extractor catalogue files will be stored
sexbgdir                 = config_data["sexbgdir"]                 # directory where Source Extractor background files will be stored

sexcheckzpdir            = config_data["sexcheckzpdir"]            # directory where Source Extractor files will be stored for zp check
sexobjcheckzpdir         = config_data["sexobjcheckzpdir"]         # directory where Source Extractor object files will be stored for zp check
sexcatcheckzpdir         = config_data["sexcatcheckzpdir"]         # directory where Source Extractor catalogue files will be stored for zp check
sexbgcheckzpdir          = config_data["sexbgcheckzpdir"]          # directory where Source Extractor background files will be stored for zp check

sexfinalcheckzpdir       = config_data["sexfinalcheckzpdir"]       # directory where Source Extractor files will be stored for final zp check
sexobjfinalcheckzpdir    = config_data["sexobjfinalcheckzpdir"]    # directory where Source Extractor object files will be stored for final zp check
sexcatfinalcheckzpdir    = config_data["sexcatfinalcheckzpdir"]    # directory where Source Extractor catalogue files will be stored for final zp check
sexbgfinalcheckzpdir     = config_data["sexbgfinalcheckzpdir"]     # directory where Source Extractor background files will be stored for final zp check

sexstackedcaldir         = config_data["sexstackedcaldir"]         # directory where Source Extractor files will be stored for stacked calibrated images
sexobjstackedcaldir      = config_data["sexobjstackedcaldir"]      # directory where Source Extractor object files will be stored for stacked calibrated images
sexcatstackedcaldir      = config_data["sexcatstackedcaldir"]      # directory where Source Extractor catalogue files will be stored for stacked calibrated images
sexbgstackedcaldir       = config_data["sexbgstackedcaldir"]       # directory where Source Extractor background files will be stored for stacked calibrated images
# log output directory
logdir                   = config_data["logdir"]                   # directory where log files will be stored

execfile(parentdir+"sourceastsex.py") # change Python virtualenv


#########################################################################################################
#                                                                                                       #
#                                       IMPORT FUNCTIONS                                                #
#                                                                                                       #
#########################################################################################################


from functools import partial  # allows extra arguments to be passed to multiprocessing (bug with lambda f'n)
import multiprocessing         # for multiprocessing
import time                    # to measure processing time
import imp                     # to force reloading

# filefuncts.py
import filefuncts
imp.reload(filefuncts)
from filefuncts import *

# calibration.py
import calibration
imp.reload(calibration)
from calibration import *

# collectzipdata.py
import collectzipdata
imp.reload(collectzipdata)
from collectzipdata import *

# callastrometry.py
import callastrometry
imp.reload(callastrometry)
from callastrometry import *

# callapass.py
#if APASS:
#    import callapass
#    imp.reload(callapass)
#    from callapass import callapass

# apassoverlap.py
#if APASS:
#    import apassoverlap
#    imp.reload(apassoverlap)
#    from apassoverlap import *

# panstarrsoverlap.py
if PANSTARRS:
    import panstarrsoverlap
    imp.reload(panstarrsoverlap)
    from panstarrsoverlap import *

# atlasoverlap.py
if ATLAS:
    import atlasoverlap
    imp.reload(atlasoverlap)
    from atlasoverlap import *

##### CHECK FOR MISSING '/' AT END OF INPUT DIRECTORIES #####

direndslash(parentdir)
direndslash(datadir)
direndslash(calfolder)
direndslash(flatfolder)
direndslash(dsffdir)
direndslash(mastercaldir)
direndslash(sciencemaskdir)
direndslash(wcsdir)
direndslash(apassdir)
direndslash(logdir)
direndslash(sexdir)
direndslash(sexobjdir)
direndslash(sexcatdir)
direndslash(sexbgdir)

if GENERATE:
    plotdir            = parentdir+"plots/"
    calplotdir         = plotdir+"calibration/"
    calibratedplotdir  = plotdir+"calibrated/"
    if not pathexists(plotdir):
        os.system("mkdir "+plotdir)
    if not pathexists(calplotdir):
        os.system("mkdir "+calplotdir)
    if not pathexists(calibratedplotdir):
        os.system("mkdir "+calibratedplotdir)


#########################################################################################################
#                                                                                                       #
#                                        OBTAIN FOLDERS                                                 #
#                                                                                                       #
#########################################################################################################


folders = []
if VERBOSE:
    print "> obtaining folders and sorting them"
for folder in os.listdir(datadir):
    folders.append(folder)
folders.sort()

badsciencefiles  = {} # dictionary of bad science images using flag input
badobjfiles      = {} # dictionary of ignored files using object input


#########################################################################################################
#                                                                                                       #
#                                          WCS SOLUTION                                                 #
#                                                                                                       #
#########################################################################################################


####################################### WCS FUNCTION ####################################################


def mainwcs(folder,folderdir,wcsfolderdir,badsciencefiles,badobjfiles,wcslogstring,file,sbool,wcslog):
    checkflag,badsciencefiles=fcheckflag(file,folder,badsciencefiles,wcslogstring,store=sbool,log=wcslog)
    if checkflag==True:
        if VERBOSE:
            print "        > skipping bad file: "+file
    else:
        checkobj,badobjfiles=fcheckobj(file,folder,badobjfiles,wcslogstring,store=sbool,log=wcslog)
        if checkobj==True:
            if VERBOSE:
                print "    > processing file "+folderdir+file
            parity = scrubwcsheader(folderdir,file,wcsfolderdir)
            callastrometry(
                wcsfolderdir+file,parity,parentdir+"astrometry.cfg",parentdir+"verify_astrometry.cfg",generate=True,filekeep=False
                )


########################################### OBTAIN WCS SOLUTION ########################################

pool_wcs = multiprocessing.Pool(ncores_wcs)

if CREATEWCS:
    if VERBOSE:
        print "\n-------------- creating world coordinate system for science images ---------------"

    # only clear out WCS directories if obtaining solution;
    # otherwise can run it once and skip this step afterward without losing files written by Astrometry.net

    if pathexists(wcsdir):
        os.system("rm -rf "+wcsdir)
    else:
        pass

    os.system("mkdir "+wcsdir)
    for folder in folders:
        os.system("mkdir "+wcsdir+folder)

    badsciencefiles  = {}
    badobjfiles      = {}
    wcslogstring     = {}

    if mode=="raw":
        sbool=False
    else:
        sbool=True

    if LOGOUTPUT:
        wcslog = True
    else:
        wcslog = False

    for folder in folders:
        print "processing folder "+folder+"..."
        folderdir = datadir+folder+"/"
        wcsfolderdir = wcsdir+folder+"/"
        files = os.listdir(folderdir)
        if len(files)>=ncores_wcs:
            mainwcs_mp = partial(mainwcs,folder,folderdir,wcsfolderdir,badsciencefiles,badobjfiles,wcslogstring,sbool=sbool,wcslog=wcslog)
            pool_wcs.map(mainwcs_mp,files)
        else:
            for file in os.listdir(folderdir):
                mainwcs(folder,folderdir,wcsfolderdir,badsciencefiles,badobjfiles,wcslogstring,file,sbool,wcslog)
else:
    if VERBOSE:
        print "\n--------------------- creating Astrometry.net WCS solution -----------------------"

#########################################################################################################
#                                                                                                       #
#                                          DATA REDUCTION                                               #
#                                                                                                       #
#########################################################################################################


############################# DATA REDUCTION FUNCTION #################################

def maincal(folder):

    start_time = time.time()

    # log output
    if LOGOUTPUT:
        logbool=True
        logstring = []
    else:
        logbool=False
        logstring = []

    missingwcs             = {} # dictionary of science images missing wcs solution
    wrongextfiles          = {} # dictionary of science images with wrong file extension
    badbiasfiles           = {} # dictionary of bad bias images using flag input
    baddarkfiles           = {} # dictionary of bad dark images using flag input
    badflatfiles           = {} # dictionary of bad flat images using flag input
    badsciencefiles        = {} # dictionary of bad science images using flag input
    badobjfiles            = {} # dictionary of ignored files using object input
    missing_flat_passband  = {} # dictionary of all files with missing flats in their respective passbands
    missing_calibration    = {} # dictionary of all files with missing calibration images
    missing_flat           = {} # dictionary of all files with missing flat images
    missing_calflat        = {} # dictionary of all files with missing calibration and flat images
    all_masterbias         = [] # list containing all master bias frames
    all_masterdarks        = {} # dictionary containing all master dark frames
    all_masterflats        = {} # dictionary containing all master flat frames
    elapsed_time           = [] # dictionary of elapsed time for each file processed

    print "------------------------- processing folder "+folder+" -------------------------"
    if LOGOUTPUT:
        logstring.append("------------------------- processing folder "+folder+" -------------------------")

    if CREATEWCS:
        folderdir = wcsdir+folder+"/"                    # folder containing wcs-corrected raw data but NOT calibration images
        if VERBOSE:
            print "calibrating folder "+folder+" using wcs-corrected raw data ..."
        if LOGOUTPUT:
            logstring.append("calibrating folder "+folder+" using wcs-corrected raw data ...")
    else:
        if pathexists(wcsdir+folder):
            folderdir = wcsdir+folder+"/"                # folder containing wcs-corrected raw data but NOT calibration images
            if VERBOSE:
                print "calibrating folder "+folder+" using wcs-corrected raw data ..."
            if LOGOUTPUT:
                logstring.append("calibrating folder "+folder+" using wcs-corrected raw data ...")
        else:
            folderdir = datadir+folder+"/"               # folder containing original raw science, bias, dark, and flat images
            if VERBOSE:
                print "calibrating folder "+folder+" using original raw data ..."
            if LOGOUTPUT:
                logstring.append("calibrating folder "+folder+" using original raw data ...")

    caldir                = datadir+folder+"/"+calfolder   # folder containing calibration images
    flatdir               = datadir+folder+"/"+flatfolder  # folder containing flat images
    dsfffolderdir         = dsffdir+folder+"/"             # folder containing dark subtracted flat fielded science images
    mastercalfolderdir    = mastercaldir+folder+"/"        # folder containing master calibration images
    sciencemaskfolderdir  = sciencemaskdir+folder+"/"      # folder containing science masks

    sbool = True

    if bothcaldirs(caldir,flatdir)==True:
        if VERBOSE:
            print "folder "+folder+" has all calibration data"
        if LOGOUTPUT:
            logstring.append("    > folder "+folder+" has all calibration data")

        # create master calibration images
        if VERBOSE:
            print "creating master bias"
        if LOGOUTPUT:
            logstring.append("    > creating master bias")
        masterbias,badbiasfiles = createmasterbias(caldir,folder,badbiasfiles,logstring,sbool,log=logbool)
        if VERBOSE:
            print "creating master darks"
        if LOGOUTPUT:
            logstring.append("    > creating master darks")
        masterdarks,baddarkfiles = createmasterdarks(caldir,masterbias,folder,baddarkfiles,logstring,sbool,log=logbool)
        if VERBOSE:
            print "creating master flats"
        if LOGOUTPUT:
            logstring.append("    > creating master flats")
        masterflats,badflatfiles = createmasterflats(flatdir,masterbias,masterdarks,folder,badflatfiles,logstring,sbool,log=logbool)

        # save master calibration images
        if VERBOSE:
            print "saving master bias"
        if LOGOUTPUT:
            logstring.append("    > saving master bias")
        savemasterbias(mastercalfolderdir,masterbias)
        if VERBOSE:
            print "saving master darks"
        if LOGOUTPUT:
            logstring.append("    > saving master darks")
        savemasterdarks(mastercalfolderdir,masterdarks)
        if VERBOSE:
            print "saving master flats"
        if LOGOUTPUT:
            logstring.append("    > saving master flats")
        savemasterflats(mastercalfolderdir,masterflats)

        # collect master calibration images
        all_masterbias.append(masterbias)
        all_masterdarks = appendmasterdarks(masterdarks,all_masterdarks)
        all_masterflats = appendmasterflats(masterflats,all_masterflats)
        sciencefiles = []
        sciencefiles = os.listdir(folderdir)
        for file in sciencefiles:
            checkobj,badobjfiles=fcheckobj(file,folder,badobjfiles,logstring,store=True,log=logbool)
            checkext,wrongextfiles=fcheckext(file,folder,".fts",wrongextfiles)
            if checkobj and checkext:
                checkflag,badsciencefiles = fcheckflag(file,folder,badsciencefiles,logstring,store=True,log=logbool)
                if checkflag:
                    if LOGOUTPUT:
                        logstring.append("    > skipping bad file: "+file)
                    if VERBOSE:
                        print "    > skipping bad file: "+file
                else:
                    if VERBOSE:
                        print "    > processing file: "+file
                    if LOGOUTPUT:
                        logstring.append("    > processing file: "+file)
                    try:
                        sciencedata, scienceheader=getdata(folderdir+file,header=True)
                        exptime=scienceheader["EXPTIME"]
                        passband=scienceheader["FILTER"]
                        if passband not in masterflats.keys():
                            if VERBOSE:
                                print file+" is missing flat in proper passband; master-master flat will be used"
                            if LOGOUTPUT:
                                logstring.append("        > "+file+" is missing flat in proper passband; master-master flat will be used")
                            missing_flat_passband = appenditemtodict(file,folder,missing_flat_passband)
                        else:
                            sciencedata_bdsff = bdsff(sciencedata,exptime,passband,masterbias,masterdarks,masterflats)
                            if sciencedata_dsff_cutoff!=None and sciencedata_bdsff.min()<sciencedata_dsff_cutoff:
                                if VERBOSE:
                                    print file+" was masked; info stored in MASK and MASKVAL header keywords"
                                if LOGOUTPUT:
                                    logstring.append("        > "+file+" was masked; info stored in MASK and MASKVAL header keywords")
                                scienceheader["MASK"]="True"
                                scienceheader["MASKVAL"]=sciencedata_dsff_cutoff
                                sciencedata_bdsff,sciencemask_bdsff = maskdata(sciencedata_bdsff,sciencedata_dsff_cutoff)
                                if GENERATE:
                                    savemask(sciencemaskfolderdir,file,sciencemask_bdsff)
                            else:
                                scienceheader["MASK"]="False"
                                scienceheader["MASKVAL"]="NA"
                            savebdsff(dsfffolderdir,file,sciencedata_bdsff,scienceheader)
                    except (IOError,OSError):
                        if VERBOSE:
                            print "        > file "+file+" missing."
                        if LOGOUTPUT:
                            logstring.append("        > file "+file+" missing (likely did not produce WCS solution).")
                        missingwcs = appenditemtodict(file,folder,missingwcs)
            else:
                if VERBOSE:
                    print "        > skipping file: "+file+" (either bad ext or missing desired object in filename)"

    elif bothcaldirs(caldir,flatdir)=="AutoFlat":
        if VERBOSE:
            print "folder "+folder+" has no flat data; master-master flats will be used"
        if LOGOUTPUT:
            logstring.append("    > folder "+folder+" has no flat data; master-master flats will be used")
        missingflat(folder,folderdir,badsciencefiles,badobjfiles,missing_flat,logstring,store=True,log=logbool)
        # create master calibration images
        if VERBOSE:
            print "creating master bias"
        if LOGOUTPUT:
            logstring.append("    > creating master bias")
        masterbias,badbiasfiles = createmasterbias(caldir,folder,badbiasfiles,logstring,sbool,log=logbool)
        if VERBOSE:
            print "creating master flats"
        if LOGOUTPUT:
            logstring.append("    > creating master darks")
        masterdarks,baddarkfiles = createmasterdarks(caldir,masterbias,folder,baddarkfiles,logstring,sbool,log=logbool)
        # save master calibration images
        if VERBOSE:
            print "saving master bias"
        if LOGOUTPUT:
            logstring.append("    > saving master bias")
        savemasterbias(mastercalfolderdir,masterbias)
        if VERBOSE:
            print "saving master darks"
        if LOGOUTPUT:
            logstring.append("    > saving master darks")
        savemasterdarks(mastercalfolderdir,masterdarks)
        # collect master calibration images
        all_masterbias.append(masterbias)
        all_masterdarks = appendmasterdarks(masterdarks,all_masterdarks)

    elif bothcaldirs(caldir,flatdir)=="Calibration":
        if VERBOSE:
            print "folder "+folder+" has no bias or dark data; master-master bias and darks will be used"
        if LOGOUTPUT:
            logstring.append("    > folder "+folder+" has no bias or dark data; master-master bias and darks will be used")
        missingcalibration(folder,folderdir,badsciencefiles,badobjfiles,missing_calibration,logstring,store=True,log=logbool)
        # create master flat images
        if VERBOSE:
            print "creating master flats"
        if LOGOUTPUT:
            logstring.append("    > creating master flats")
        masterflats,badflatfiles = createmasterflats(flatdir,masterbias,masterdarks,folder,badflatfiles,logstring,sbool,log=logbool)
        # save master calibration images
        if VERBOSE:
            print "saving master flats"
        if LOGOUTPUT:
            logstring.append("    > saving master flats")
        savemasterflats(mastercalfolderdir,masterflats)
        # collect master flat images
        all_masterflats = appendmasterflats(masterflats,all_masterflats)

    elif bothcaldirs(caldir,flatdir)==False:
        if LOGOUTPUT:
            logstring.append("no bias, dark, or flat data for folder "+folder+"; master-master images will be used")
        if VERBOSE:
            print "no bias, dark, or flat data for folder "+folder+"; master-master images will be used"
        missingcalflat(folder,folderdir,badsciencefiles,badobjfiles,missing_calflat,logstring,store=True,log=logbool)

    end_time = time.time()
    elapsed_time_i = end_time - start_time
    elapsed_time.append(elapsed_time_i)

    if LOGOUTPUT:
        logstring.append("elapsed time: "+str(elapsed_time)+" seconds\n")

    return badbiasfiles,baddarkfiles,badflatfiles,badsciencefiles,missing_flat_passband,missing_calibration,missing_flat,missing_calflat,all_masterbias,all_masterdarks,all_masterflats,wrongextfiles,elapsed_time,logstring


############################## PERFORM DATA REDUCTION  ##############################

if sciencedata_dsff_cutoff!=None:
    sciencedata_dsff_cutoff = float(sciencedata_dsff_cutoff)

checkdirectories = [mastercaldir,MMcaldir,dsffdir,sciencemaskdir,logdir]

pool_cal = multiprocessing.Pool(ncores_cal)

if REDUCE:
    if VERBOSE:
        print "\n------------------------- proceeding with data reduction -------------------------"

    # only clear out directories if performing calibration
    for dir in checkdirectories:
        if pathexists(dir):
            os.system("rm -rf "+dir)
        else:
            pass
        os.system("mkdir "+dir)
        if VERBOSE:
            print "cleared directory: "+dir
        if dir!= MMcaldir and dir!=logdir:
            for folder in folders:
                os.system("mkdir "+dir+folder)
                if VERBOSE:
                    print "cleared directory: "+dir+folder

    badbiasfiles_zip,baddarkfiles_zip,badflatfiles_zip,badsciencefiles_zip,missing_flat_passband_zip,missing_calibration_zip,missing_flat_zip,missing_calflat_zip,all_masterbias_zip,all_masterdarks_zip,all_masterflats_zip,wrongextfiles_zip,elapsed_time_zip,logstring_zip = zip(*pool_cal.map(maincal,folders))

    # fix all zipped lists and dictionaries to regular lists and dictionaries

    badbiasfiles           = collectzipdicts(badbiasfiles_zip)
    baddarkfiles           = collectzipdicts(baddarkfiles_zip)
    badflatfiles           = collectzipdicts(badflatfiles_zip)
    badsciencefiles        = collectzipdicts(badsciencefiles_zip)
    missing_flat_passband  = collectzipdicts(missing_flat_passband_zip)
    missing_calibration    = collectzipdicts(missing_calibration_zip)
    missing_flat           = collectzipdicts(missing_flat_zip)
    missing_calflat        = collectzipdicts(missing_calflat_zip)
    wrongextfiles          = collectzipdicts(wrongextfiles_zip)
    elapsed_time           = collectziplists(elapsed_time_zip)
    logstring              = collectziplists(logstring_zip)

    if VERBOSE:
        print "\n--------------- creating master-master bias, dark, and flat images ---------------"
    if LOGOUTPUT:
        logstring.append("\n--------------- creating master-master bias, dark, and flat images ---------------")

    if VERBOSE:
        print "creating master-master bias"
    if LOGOUTPUT:
        logstring.append("creating master-master bias")
    all_masterbias = collectziplists(all_masterbias_zip)
    MMbias = createMMbias(all_masterbias)

    if VERBOSE:
        print "creating master-master darks"
    if LOGOUTPUT:
        logstring.append("creating master-master darks")
    all_masterdarks = collectzipdicts(all_masterdarks_zip)
    MMdarks = createMMdarks(all_masterdarks)

    if VERBOSE:
        print "creating master-master flats"
    if LOGOUTPUT:
        logstring.append("creating master-master flats")
    all_masterflats = collectzipdicts(all_masterflats_zip)
    MMflats = createMMflats(all_masterflats)

    # save master-master images
    if VERBOSE:
        print "saving master-master bias"
    if LOGOUTPUT:
        logstring.append("saving master-master bias")
    savemasterbias(MMcaldir,MMbias)
    if VERBOSE:
        print "saving master-master darks"
    if LOGOUTPUT:
        logstring.append("saving master-master darks")
    savemasterdarks(MMcaldir,MMdarks)
    if VERBOSE:
        print "saving master-master flats"
    if LOGOUTPUT:
        logstring.append("saving master-master flats")
    savemasterflats(MMcaldir,MMflats)

    if VERBOSE:
        print "\n------------------ checking for missing biases, darks, and flats -----------------"
    if LOGOUTPUT:
        logstring.append("\n------------------ checking for missing biases, darks, and flats -----------------")

    if LOGOUTPUT:
        logbool=True
    else:
        logbool=False

    # missing flats
    if not missing_flat:
        if VERBOSE:
            print "there are no science images missing flat frames"
        if LOGOUTPUT:
            logstring.append("there are no science images missing flat frames")
    else:
        if VERBOSE:
            print "there are science images missing all flat frames"
        if LOGOUTPUT:
            logstring.append("there are science images missing all flat frames")
        badsciencefiles,badobjfiles,wrongextfiles,logstring=bdsff_missingdata(
            missing_flat,VERBOSE,GENERATE,logstring,badsciencefiles,badobjfiles,wrongextfiles,log=logbool,useMMflat=True
            )

    # missing flats in proper passband
    if not missing_flat_passband:
        if VERBOSE:
            print "there are no science images missing flats in their proper passbands"
        if LOGOUTPUT:
            logstring.append("there are no science images missing flats in their proper passbands")
    else:
        if VERBOSE:
            print "there are science images missing flat frames in their proper passbands"
        if LOGOUTPUT:
            logstring.append("there are science images missing flat frames in their proper passbands")
        badsciencefiles,badobjfiles,wrongextfiles,logstring=bdsff_missingdata(
            missing_flat_passband,VERBOSE,GENERATE,logstring,badsciencefiles,badobjfiles,wrongextfiles,log=logbool,useMMflat=True
            )

    # missing bias, darks
    if not missing_calibration:
        if VERBOSE:
            print "there are no science images missing all bias and dark frames"
        if LOGOUTPUT:
            logstring.append("there are no science images missing all bias and dark frames")
    else:
        if VERBOSE:
            print "there are science images missing all bias and dark frames"
        if LOGOUTPUT:
            logstring.append("there are science images missing all bias and dark frames")
        bdsff_missingdata(missing_calibration,VERBOSE,GENERATE,logstring,log=logbool,useMMbias=True,useMMdark=True)

    # missing flat + bias and darks
    if not missing_calflat:
        if VERBOSE:
            print "there are no science images missing all bias, dark, and flat frames"
        if LOGOUTPUT:
            logstring.append("there are no science images missing all bias, dark, and flat frames")
    else:
        if VERBOSE:
            print "there are science images missing all bias, dark, and flat frames"
        if LOGOUTPUT:
            logstring.append("there are science images missing all bias, dark, and flat frames")
        badsciencefiles,badobjfiles,wrongextfiles,logstring=bdsff_missingdata(
            missing_calflat,VERBOSE,GENERATE,logstring,badsciencefiles,badobjfiles,wrongextfiles,log=logbool,useMMbias=True,useMMdark=True,useMMflat=True
            )

    if LOGOUTPUT:

        # write calibration log
        if VERBOSE:
            print "writing calibration log"
        f=open(logdir+"calibration.log","w+")
        for item in logstring:
            f.write("%s\n" % item)
        f.close()

        # write bad files log
        if VERBOSE:
            print "> writing bad file log"
        f=open(logdir+"badfiles.log","w+")
        filedicts = [badbiasfiles,baddarkfiles,badflatfiles,badsciencefiles,badobjfiles]
        dictstrings = ["bias files flagged as bad","dark files flagged as being bad","flat files flagged as being bad","science files flagged as being bad","science files lacking desired object in filename"]
        for folder in folders:
            f.write("%s\n" % folder)
            for i in range(len(filedicts)):
                filedict = filedicts[i]
                dictstring = dictstrings[i]
                f.write("%s\n" % dictstring)
                if folder in filedict.keys():
                    filelist = filedict[folder]
                    if len(filelist)==0:
                        f.write("None\n")
                    elif len(filelist)!=0:
                        for file in filelist:
                            f.write("%s\n" % file)
                else:
                    f.write("None\n")
            f.write("\n")
    f.close()

else:
    if VERBOSE:
        print "\n----------------------------- skipping data reduction ----------------------------"

#########################################################################################################
#                                                                                                       #
#                                     PAN-STARRS CALIBRATION                                            #
#                                                                                                       #
#########################################################################################################

#########################################################################################################
#                                                                                                       #
#                                     ZEROPOINT CALIBRATION                                             #
#                                                                                                       #
#########################################################################################################

######################################## ZEROPOINT FUNCTION #############################################

panstarrspassbands = ["g","r","i","z","y"]

def mainzeropointcal_panstarrs(panstarrsdata,folder,folderdir,file,logbool=False,CHECKZP=False):
    #
    # Measures the magnitude zero point by comparing Source Extractor measurements to Pan-STARRS. Uses the slope of
    # (Pan-STARRS mag - SExtractor mag) vs. Pan-STARRS mag to determine if the image will be included; a non-zero slope
    # within 3-sigma uncertainties will be returned in the 'badcalfiles' dictionary and not used any further.
    #

    logstring                              = []

    ######################## long exposure (120 s) ########################
    # zero slopes and 'good' zero point values
    stackscienceimages_longexp_dict        = {}
    zeropoints_longexp_dict                = {}
    zeropoint_errs_longexp_dict            = {}
    airmasses_longexp_dict                 = {}
    timeobs_longexp_dict                   = {}
    FWHM_longexp_dict                      = {}
    # matched and cleaned
    mag_sex_matches_clean_longexp_dict     = {}
    magerr_sex_matches_clean_longexp_dict  = {}
    mag_ps_matches_clean_longexp_dict      = {}
    magerr_ps_matches_clean_longexp_dict   = {}
    # matched
    mag_sex_matches_longexp_dict           = {}
    magerr_sex_matches_longexp_dict        = {}
    mag_ps_matches_longexp_dict            = {}
    magerr_ps_matches_longexp_dict         = {}
    # all zero slopes ('good' and 'bad' zero point values) for plotting
    zeropoints_all_longexp_dict            = {}
    zeropoints_err_all_longexp_dict        = {}
    airmasses_all_longexp_dict             = {}
    timeobs_all_longexp_dict               = {}
    FWHM_all_longexp_dict                  = {}
    # zero slopes but 'poor' zeropoint values AND non-zero slopes
    donotstackscienceimages_longexp_dict   = {}

    ######################## short exposure (5 s) ########################
    # zero slopes and 'good' zero point values
    stackscienceimages_shortexp_dict        = {}
    zeropoints_shortexp_dict                = {}
    zeropoint_errs_shortexp_dict            = {}
    airmasses_shortexp_dict                 = {}
    timeobs_shortexp_dict                   = {}
    FWHM_shortexp_dict                      = {}
    # matched and cleaned
    mag_sex_matches_clean_shortexp_dict     = {}
    magerr_sex_matches_clean_shortexp_dict  = {}
    mag_ps_matches_clean_shortexp_dict      = {}
    magerr_ps_matches_clean_shortexp_dict   = {}
    # matched
    mag_sex_matches_shortexp_dict           = {}
    magerr_sex_matches_shortexp_dict        = {}
    mag_ps_matches_shortexp_dict            = {}
    magerr_ps_matches_shortexp_dict         = {}
    # all zero slopes ('good' and 'bad' zero point values) for plotting
    zeropoints_all_shortexp_dict            = {}
    zeropoints_err_all_shortexp_dict        = {}
    airmasses_all_shortexp_dict             = {}
    timeobs_all_shortexp_dict               = {}
    FWHM_all_shortexp_dict                  = {}
    # zero slopes but 'poor' zeropoint values AND non-zero slopes
    donotstackscienceimages_shortexp_dict   = {}

    sciencedata   = fits.getdata(folderdir+file)
    scienceheader = fits.getheader(folderdir+file)
    passband      = scienceheader["FILTER"]
    airmass       = scienceheader["AIRMASS"]
    exptime       = scienceheader["EXPTIME"]
    timeobs       = scienceheader["TIME-OBS"]
    timeobs_24h   = ftimeobs_to24h(timeobs)
    FWHM          = scienceheader["FWHM"]
    object        = re.split("-|\.",file)[0]
    if passband in panstarrspassbands:
        if VERBOSE:
            print "    > processing file: "+file
        if LOGOUTPUT:
            logstring.append("    > processing file: "+str(file))

        # run Source Extractor if indicated otherwise open already existing file
        if CHECKZP==False:
            sexcat = sexcall(file,folder,folderdir,sexobjdir,sexcatdir,sexbgdir,sexbool=SEXTRACTOR,CHECKZP=CHECKZP)
        elif CHECKZP==True:
            sexcat = sexcall(file,folder,folderdir,sexobjcheckzpdir,sexcatcheckzpdir,sexbgcheckzpdir,sexbool=SEXTRACTOR,CHECKZP=CHECKZP)

        # find Pan-STARRS overlap with SExtractor
        ZEROSLOPE,zeropoint,zeropoint_error,GOODZP,CHECKOBSTIME,mag_sex_matches_clean,magerr_sex_matches_clean,mag_ps_matches_clean,magerr_ps_matches_clean,mag_sex_matches,magerr_sex_matches,mag_ps_matches,magerr_ps_matches = magcalibration(file,folder,folderdir,panstarrsdata,sexcat,VERBOSE,GENERATE,CHECKZP)

        ######################## long exposure (120 s) ########################
        if exptime==120.:
            print passband, exptime, "(long exposure)"
            # collect files with zero slopes and 'good' zero point values
            if (ZEROSLOPE==True) and (GOODZP==True) and (CHECKOBSTIME==True):
                # zero slopes and 'good' zero point values
                stackscienceimages_longexp_dict          = appenditemtodict(folderdir+file,object,stackscienceimages_longexp_dict)
                zeropoints_longexp_dict                  = appenditemtodict(zeropoint,passband,zeropoints_longexp_dict)
                zeropoint_errs_longexp_dict              = appenditemtodict(zeropoint_error,passband,zeropoint_errs_longexp_dict)
                airmasses_longexp_dict                   = appenditemtodict(airmass,passband,airmasses_longexp_dict)
                timeobs_longexp_dict                     = appenditemtodict(timeobs_24h,passband,timeobs_longexp_dict)
                FWHM_longexp_dict                        = appenditemtodict(FWHM,passband,FWHM_longexp_dict)
                # matched and cleaned
                mag_sex_matches_clean_longexp_dict       = appenditemtodict(mag_sex_matches_clean,passband,mag_sex_matches_clean_longexp_dict)
                magerr_sex_matches_clean_longexp_dict    = appenditemtodict(magerr_sex_matches_clean,passband,magerr_sex_matches_clean_longexp_dict)
                mag_ps_matches_clean_longexp_dict        = appenditemtodict(mag_ps_matches_clean,passband,mag_ps_matches_clean_longexp_dict)
                magerr_ps_matches_clean_longexp_dict     = appenditemtodict(magerr_ps_matches_clean,passband,magerr_ps_matches_clean_longexp_dict)
                # matched
                mag_sex_matches_longexp_dict             = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_longexp_dict)
                magerr_sex_matches_longexp_dict          = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_longexp_dict)
                mag_ps_matches_longexp_dict              = appenditemtodict(mag_ps_matches,passband,mag_ps_matches_longexp_dict)
                magerr_ps_matches_longexp_dict           = appenditemtodict(magerr_ps_matches,passband,magerr_ps_matches_longexp_dict)
                if CHECKZP==False:
                    # add zero-point to FITS header
                    fits.setval(folderdir+file, "ZP", value=zeropoint)
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=True)
                    # write calibrated image
                    calibratedfolderdir    = pscalibrateddir+folder+"/"
                    newfilename            = file.split(".fts")[0]+"_calibrated.fts"
                    #sciencedata_cal       = fZPcorrADUs(sciencedata,zeropoint)
                    if VERBOSE:
                        print "writing calibrated image: "+calibratedfolderdir+newfilename
                    fits.writeto(calibratedfolderdir+newfilename,sciencedata,scienceheader,overwrite=True)
            # collect files with non-zero slopes and 'bad' zero-point values
            else:
                donotstackscienceimages_longexp_dict     = appenditemtodict(folderdir+file,object,donotstackscienceimages_longexp_dict)
                if zeropoint!=None:
                    zeropoints_all_longexp_dict          = appenditemtodict(zeropoint,passband,zeropoints_all_longexp_dict)
                    zeropoints_err_all_longexp_dict      = appenditemtodict(zeropoint_error,passband,zeropoints_err_all_longexp_dict)
                    airmasses_all_longexp_dict           = appenditemtodict(airmass,passband,airmasses_all_longexp_dict)
                    timeobs_all_longexp_dict             = appenditemtodict(timeobs_24h,passband,timeobs_all_longexp_dict)
                    FWHM_all_longexp_dict                = appenditemtodict(FWHM,passband,FWHM_all_longexp_dict)
                    # matched
                    mag_sex_matches_longexp_dict         = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_longexp_dict)
                    magerr_sex_matches_longexp_dict      = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_longexp_dict)
                    mag_ps_matches_longexp_dict          = appenditemtodict(mag_ps_matches,passband,mag_ps_matches_longexp_dict)
                    magerr_ps_matches_longexp_dict       = appenditemtodict(magerr_ps_matches,passband,magerr_ps_matches_longexp_dict)
                if CHECKZP==False:
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=False)
        
        ######################## short exposure (5 s) ########################
        elif exptime==5.:
            # collect files with zero slopes and 'good' zero point values
            if (ZEROSLOPE==True) and (GOODZP==True) and (CHECKOBSTIME==True):
                # zero slopes and 'good' zero point values
                stackscienceimages_shortexp_dict          = appenditemtodict(folderdir+file,object,stackscienceimages_shortexp_dict)
                zeropoints_shortexp_dict                  = appenditemtodict(zeropoint,passband,zeropoints_shortexp_dict)
                zeropoint_errs_shortexp_dict              = appenditemtodict(zeropoint_error,passband,zeropoint_errs_shortexp_dict)
                airmasses_shortexp_dict                   = appenditemtodict(airmass,passband,airmasses_shortexp_dict)
                timeobs_shortexp_dict                     = appenditemtodict(timeobs_24h,passband,timeobs_shortexp_dict)
                FWHM_shortexp_dict                        = appenditemtodict(FWHM,passband,FWHM_shortexp_dict)
                # matched and cleaned
                mag_sex_matches_clean_shortexp_dict       = appenditemtodict(mag_sex_matches_clean,passband,mag_sex_matches_clean_shortexp_dict)
                magerr_sex_matches_clean_shortexp_dict    = appenditemtodict(magerr_sex_matches_clean,passband,magerr_sex_matches_clean_shortexp_dict)
                mag_ps_matches_clean_shortexp_dict        = appenditemtodict(mag_ps_matches_clean,passband,mag_ps_matches_clean_shortexp_dict)
                magerr_ps_matches_clean_shortexp_dict     = appenditemtodict(magerr_ps_matches_clean,passband,magerr_ps_matches_clean_shortexp_dict)
                # matched
                mag_sex_matches_shortexp_dict             = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_shortexp_dict)
                magerr_sex_matches_shortexp_dict          = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_shortexp_dict)
                mag_ps_matches_shortexp_dict              = appenditemtodict(mag_ps_matches,passband,mag_ps_matches_shortexp_dict)
                magerr_ps_matches_shortexp_dict           = appenditemtodict(magerr_ps_matches,passband,magerr_ps_matches_shortexp_dict)
                if CHECKZP==False:
                    # add zero-point to FITS header
                    fits.setval(folderdir+file, "ZP", value=zeropoint)
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=True)
                    # write calibrated image
                    calibratedfolderdir = pscalibrateddir+folder+"/"
                    newfilename         = file.split(".fts")[0]+"_calibrated.fts"
                    #sciencedata_cal    = fZPcorrADUs(sciencedata,zeropoint)
                    if VERBOSE:
                        print "writing calibrated image: "+calibratedfolderdir+newfilename
                    fits.writeto(calibratedfolderdir+newfilename,sciencedata,scienceheader,overwrite=True)
            # collect files with non-zero slopes and 'bad' zero-point values
            else:
                donotstackscienceimages_shortexp_dict     = appenditemtodict(folderdir+file,object,donotstackscienceimages_shortexp_dict)
                # matched
                mag_sex_matches_shortexp_dict             = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_shortexp_dict)
                magerr_sex_matches_shortexp_dict          = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_shortexp_dict)
                mag_ps_matches_shortexp_dict              = appenditemtodict(mag_ps_matches,passband,mag_ps_matches_shortexp_dict)
                magerr_ps_matches_shortexp_dict           = appenditemtodict(magerr_ps_matches,passband,magerr_ps_matches_shortexp_dict)
                if zeropoint!=None:
                    # bad zero-points
                    zeropoints_all_shortexp_dict          = appenditemtodict(zeropoint,passband,zeropoints_all_shortexp_dict)
                    zeropoints_err_all_shortexp_dict      = appenditemtodict(zeropoint_error,passband,zeropoints_err_all_shortexp_dict)
                    airmasses_all_shortexp_dict           = appenditemtodict(airmass,passband,airmasses_all_shortexp_dict)
                    timeobs_all_shortexp_dict             = appenditemtodict(timeobs_24h,passband,timeobs_all_shortexp_dict)
                    FWHM_all_shortexp_dict                = appenditemtodict(FWHM,passband,FWHM_all_shortexp_dict)
                    # matched
                    mag_sex_matches_shortexp_dict         = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_shortexp_dict)
                    magerr_sex_matches_shortexp_dict      = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_shortexp_dict)
                    mag_ps_matches_shortexp_dict          = appenditemtodict(mag_ps_matches,passband,mag_ps_matches_shortexp_dict)
                    magerr_ps_matches_shortexp_dict       = appenditemtodict(magerr_ps_matches,passband,magerr_ps_matches_shortexp_dict)
                if CHECKZP==False:
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=False)

    return logstring,stackscienceimages_longexp_dict,zeropoints_longexp_dict,zeropoint_errs_longexp_dict,airmasses_longexp_dict,timeobs_longexp_dict,FWHM_longexp_dict,mag_sex_matches_clean_longexp_dict,magerr_sex_matches_clean_longexp_dict,mag_ps_matches_clean_longexp_dict,magerr_ps_matches_clean_longexp_dict,mag_sex_matches_longexp_dict,magerr_sex_matches_longexp_dict,mag_ps_matches_longexp_dict,magerr_ps_matches_longexp_dict,zeropoints_all_longexp_dict,zeropoints_err_all_longexp_dict,airmasses_all_longexp_dict,timeobs_all_longexp_dict,FWHM_all_longexp_dict,donotstackscienceimages_longexp_dict,stackscienceimages_shortexp_dict,zeropoints_shortexp_dict,zeropoint_errs_shortexp_dict,airmasses_shortexp_dict,timeobs_shortexp_dict,FWHM_shortexp_dict,mag_sex_matches_clean_shortexp_dict,magerr_sex_matches_clean_shortexp_dict,mag_ps_matches_clean_shortexp_dict,magerr_ps_matches_clean_shortexp_dict,mag_sex_matches_shortexp_dict,magerr_sex_matches_shortexp_dict,mag_ps_matches_shortexp_dict,magerr_ps_matches_shortexp_dict,zeropoints_all_shortexp_dict,zeropoints_err_all_shortexp_dict,airmasses_all_shortexp_dict,timeobs_all_shortexp_dict,FWHM_all_shortexp_dict,donotstackscienceimages_shortexp_dict

if CALIBRATE and PANSTARRS:
    CHECKZP=False

    # read in catalogue data
    panstarrscatf = panstarrsdir+panstarrscat

    if VERBOSE:
        print "\n---------------------- proceeding with Pan-STARRS calilbration -------------------"

    print "Reading in Pan-STARRS catalogue..."
    hdul_panstarrs       = fits.open(panstarrscatf)
    panstarrsdata        = hdul_panstarrs[1].data
    PS_objInfoFlag       = panstarrsdata["objInfoFlag"]       # Information flag bitmask indicating details of the photometry. Values listed in ObjectInfoFlags.
    PS_qualityFlag       = panstarrsdata["qualityFlag"]       # Subset of objInfoFlag denoting whether this object is real or a likely false positive. Values listed in ObjectQualityFlags.
    PS_raMean            = panstarrsdata["raMean"]            # Right ascension from single epoch detections (weighted mean) in equinox J2000 at the mean epoch given by epochMean.
    PS_decMean           = panstarrsdata["decMean"]           # Declination from single epoch detections (weighted mean) in equinox J2000 at the mean epoch given by epochMean.
    PS_raMeanErr         = panstarrsdata["raMeanErr"]         # Right ascension standard deviation from single epoch detections.
    PS_decMeanErr        = panstarrsdata["decMeanErr"]        # Declination standard deviation from single epoch detections.
    PS_gMeanPSFMag       = panstarrsdata["gMeanPSFMag"]       # Mean PSF magnitude from g filter detections.
    PS_gMeanPSFMagErr    = panstarrsdata["gMeanPSFMagErr"]    # Error in mean PSF magnitude from g filter detections.
    PS_gMeanKronMag      = panstarrsdata["gMeanKronMag"]      # Mean Kron (1980) magnitude from g filter detections.
    PS_gMeanKronMagErr   = panstarrsdata["gMeanKronMagErr"]   # Error in mean Kron (1980) magnitude from g filter detections.
    PS_gMeanApMag        = panstarrsdata["gMeanApMag"]        # Mean aperture magnitude from g filter detections.
    PS_gMeanApMagErr     = panstarrsdata["gMeanApMagErr"]     # Error in mean aperture magnitude from g filter detections.
    PS_gFlags            = panstarrsdata["gFlags"]            # Information flag bitmask for mean object from g filter detections. Values listed in ObjectFilterFlags.
    PS_rQfPerfect        = panstarrsdata["rQfPerfect"]        # Maximum PSF weighted fraction of pixels totally unmasked from r filter detections.
    PS_rMeanPSFMag       = panstarrsdata["rMeanPSFMag"]       # Mean PSF magnitude from r filter detections.
    PS_rMeanPSFMagErr    = panstarrsdata["rMeanPSFMagErr"]    # Error in mean PSF magnitude from r filter detections.
    PS_rMeanKronMag      = panstarrsdata["rMeanKronMag"]      # Mean Kron (1980) magnitude from r filter detections.
    PS_rMeanKronMagErr   = panstarrsdata["rMeanKronMagErr"]   # Error in mean Kron (1980) magnitude from r filter detections.
    PS_rMeanApMag        = panstarrsdata["rMeanApMag"]        # Mean aperture magnitude from r filter detections.
    PS_rMeanApMagErr     = panstarrsdata["rMeanApMagErr"]     # Error in mean aperture magnitude from r filter detections.
    PS_rFlags            = panstarrsdata["rFlags"]            # Information flag bitmask for mean object from r filter detections. Values listed in ObjectFilterFlags.
    PS_iQfPerfect        = panstarrsdata["iQfPerfect"]        # Maximum PSF weighted fraction of pixels totally unmasked from i filter detections.
    PS_iMeanPSFMag       = panstarrsdata["iMeanPSFMag"]       # Mean PSF magnitude from i filter detections.
    PS_iMeanPSFMagErr    = panstarrsdata["iMeanPSFMagErr"]    # Error in mean PSF magnitude from i filter detections.
    PS_iMeanKronMag      = panstarrsdata["iMeanKronMag"]      # Mean Kron (1980) magnitude from i filter detections.
    PS_iMeanKronMagErr   = panstarrsdata["iMeanKronMagErr"]   # Error in mean Kron (1980) magnitude from i filter detections.
    PS_iMeanApMag        = panstarrsdata["iMeanApMag"]        # Mean aperture magnitude from i filter detections.
    PS_iMeanApMagErr     = panstarrsdata["iMeanApMagErr"]     # Error in mean aperture magnitude from i filter detections.
    PS_iFlags            = panstarrsdata["iFlags"]            # Information flag bitmask for mean object from i filter detections. Values listed in ObjectFilterFlags.
    PS_zQfPerfect        = panstarrsdata["zQfPerfect"]        # Maximum PSF weighted fraction of pixels totally unmasked from z filter detections.
    PS_zMeanPSFMag       = panstarrsdata["zMeanPSFMag"]       # Mean PSF magnitude from z filter detections.
    PS_zMeanPSFMagErr    = panstarrsdata["zMeanPSFMagErr"]    # Error in mean PSF magnitude from z filter detections.
    PS_zMeanKronMag      = panstarrsdata["zMeanKronMag"]      # Mean Kron (1980) magnitude from z filter detections.
    PS_zMeanKronMagErr   = panstarrsdata["zMeanKronMagErr"]   # Error in mean Kron (1980) magnitude from z filter detections.
    PS_zMeanApMag        = panstarrsdata["zMeanApMag"]        # Mean aperture magnitude from z filter detections.
    PS_zMeanApMagErr     = panstarrsdata["zMeanApMagErr"]     # Error in mean aperture magnitude from z filter detections.
    PS_zFlags            = panstarrsdata["zFlags"]            # Information flag bitmask for mean object from z filter detections. Values listed in ObjectFilterFlags.
    PS_yQfPerfect        = panstarrsdata["yQfPerfect"]        # Maximum PSF weighted fraction of pixels totally unmasked from y filter detections.
    PS_yMeanPSFMag       = panstarrsdata["yMeanPSFMag"]       # Mean PSF magnitude from y filter detections.
    PS_yMeanPSFMagErr    = panstarrsdata["yMeanPSFMagErr"]    # Error in mean PSF magnitude from y filter detections.
    PS_yMeanKronMag      = panstarrsdata["yMeanKronMag"]      # Mean Kron (1980) magnitude from y filter detections.
    PS_yMeanKronMagErr   = panstarrsdata["yMeanKronMagErr"]   # Error in mean Kron (1980) magnitude from y filter detections.
    PS_yMeanApMag        = panstarrsdata["yMeanApMag"]        # Mean aperture magnitude from y filter detections.
    PS_yMeanApMagErr     = panstarrsdata["yMeanApMagErr"]     # Error in mean aperture magnitude from y filter detections.
    PS_yFlags            = panstarrsdata["yFlags"]            # Information flag bitmask for mean object from y filter detections. Values listed in ObjectFilterFlags.
    # APERTURE PHOTOMETRY
    if panstarrsphot=="Ap":
        PS_rMag        = panstarrsdata["rMeanApMag"]
        PS_rMagErr     = panstarrsdata["rMeanApMagErr"]
        PS_iMag        = panstarrsdata["iMeanApMag"]
        PS_iMagErr     = panstarrsdata["iMeanApMagErr"]
        PS_zMag        = panstarrsdata["zMeanApMag"]
        PS_zMagErr     = panstarrsdata["zMeanApMagErr"]
    # PSF PHOTOMETRY
    elif panstarrsphot=="PSF":
        PS_rMag       = panstarrsdata["rMeanPSFMag"]
        PS_rMagErr    = panstarrsdata["rMeanPSFMagErr"]
        PS_iMag       = panstarrsdata["iMeanPSFMag"]
        PS_iMagErr    = panstarrsdata["iMeanPSFMagErr"]
        PS_zMag       = panstarrsdata["zMeanPSFMag"]
        PS_zMagErr    = panstarrsdata["zMeanPSFMagErr"]
    # KRON PHOTOMETRY
    elif panstarrsphot=="Kron":
        PS_rMag      = panstarrsdata["rMeanKronMag"]
        PS_rMagErr   = panstarrsdata["rMeanKronMagErr"]
        PS_iMag      = panstarrsdata["iMeanKronMag"]
        PS_iMagErr   = panstarrsdata["iMeanKronMagErr"]
        PS_zMag      = panstarrsdata["zMeanKronMag"]
        PS_zMagErr   = panstarrsdata["zMeanKronMagErr"]

    # trimmed version of panstarrs data
    panstarrstrimcatf = panstarrsdir+panstarrstrimcat
    t = Table([PS_raMean,PS_raMeanErr,PS_decMean,PS_decMeanErr,PS_rMag,PS_rMagErr,PS_rFlags,PS_iMag,PS_iMagErr,PS_iFlags,PS_zMag,PS_zMagErr,PS_zFlags,PS_objInfoFlag,PS_qualityFlag],names=("raMean","raMeanErr","decMean","decMeanErr","rMean"+panstarrsphot+"Mag","rMean"+panstarrsphot+"MagErr","rFlags","iMean"+panstarrsphot+"Mag","iMean"+panstarrsphot+"MagErr","iFlags","zMean"+panstarrsphot+"Mag","zMean"+panstarrsphot+"MagErr","zFlags","objInfoFlag","qualityFlag"))
    t.write(panstarrstrimcatf,format="fits",overwrite=True)

    hdul_panstarrstrim       = fits.open(panstarrstrimcatf)
    panstarrsdata_trim       = hdul_panstarrstrim[1].data

    #########################################################################################################
    #                                                                                                       #
    #                                              NO NULLS                                                 #
    #                                            NO SATURATION                                              #
    #                                              NO PAIRS                                                 #
    #                                                                                                       #
    #########################################################################################################

    if VERBOSE:
        "Applying quality cut on Pan-STARRS data..."

    # r band
    ii                      = (PS_rMag!=-999.0) & (PS_rMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_rFlags<8388608.)
    PS_rMag_qual            = PS_rMag[ii]
    PS_rMagErr_qual         = PS_rMagErr[ii]
    # i band
    ii                      = (PS_iMag!=-999.0) & (PS_iMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_rFlags<8388608.)
    PS_iMag_qual            = PS_iMag[ii]
    PS_iMagErr_qual         = PS_iMagErr[ii]
    # z band
    ii                      = (PS_zMag!=-999.0) & (PS_zMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_rFlags<8388608.)
    PS_zMag_qual            = PS_zMag[ii]
    PS_zMagErr_qual         = PS_zMagErr[ii]

    # r band: PSF
    ii                      = (PS_rMeanPSFMag!=-999.0) & (PS_rMeanPSFMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_rFlags<8388608.)
    PS_rMeanPSFMag_qual     = PS_rMeanPSFMag[ii]
    PS_rMeanPSFMagErr_qual  = PS_rMeanPSFMagErr[ii]
    # r band: Kron
    ii                      = (PS_rMeanKronMag!=-999.0) & (PS_rMeanKronMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_rFlags<8388608.)
    PS_rMeanKronMag_qual    = PS_rMeanKronMag[ii]
    PS_rMeanKronMagErr_qual = PS_rMeanKronMagErr[ii]
    # r band: Aperture
    ii                      = (PS_rMeanApMag!=-999.0) & (PS_rMeanApMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_rFlags<8388608.)
    PS_r_qual               = PS_rMeanApMag[ii]
    PS_rerr_qual            = PS_rMeanApMagErr[ii]

    # i band: PSF
    ii                      = (PS_iMeanPSFMag!=-999.0) & (PS_iMeanPSFMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_iFlags<8388608.)
    PS_iMeanPSFMag_qual     = PS_iMeanPSFMag[ii]
    PS_iMeanPSFMagErr_qual  = PS_iMeanPSFMagErr[ii]
    # i band: Kron
    ii                      = (PS_iMeanKronMag!=-999.0) & (PS_iMeanKronMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_iFlags<8388608.)
    PS_iMeanKronMag_qual    = PS_iMeanKronMag[ii]
    PS_iMeanKronMagErr_qual = PS_iMeanKronMagErr[ii]
    # i band: Aperture
    ii                      = (PS_iMeanApMag!=-999.0) & (PS_iMeanApMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_iFlags<8388608.)
    PS_i_qual               = PS_iMeanApMag[ii]
    PS_ierr_qual            = PS_iMeanApMagErr[ii]

    # z band: PSF
    ii                      = (PS_zMeanPSFMag!=-999.0) & (PS_zMeanPSFMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_zFlags<8388608.)
    PS_zMeanPSFMag_qual     = PS_zMeanPSFMag[ii]
    PS_zMeanPSFMagErr_qual  = PS_zMeanPSFMagErr[ii]
    # z band: Kron
    ii                      = (PS_zMeanKronMag!=-999.0) & (PS_zMeanKronMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_zFlags<8388608.)
    PS_zMeanKronMag_qual    = PS_zMeanKronMag[ii]
    PS_zMeanKronMagErr_qual = PS_zMeanKronMagErr[ii]
    # z band: Aperture
    ii                      = (PS_zMeanApMag!=-999.0) & (PS_rMeanApMagErr!=-999.0) & (PS_qualityFlag<128.) & (PS_zFlags<8388608.)
    PS_z_qual               = PS_zMeanApMag[ii]
    PS_zerr_qual            = PS_zMeanApMagErr[ii]

    #                              PSF                                        Kron                                         Aperture
    panstarrsdata_mags = np.array([PS_rMeanPSFMag_qual,PS_rMeanPSFMagErr_qual,PS_rMeanKronMag_qual,PS_rMeanKronMagErr_qual,PS_rMeanApMag,PS_rMeanApMagErr,
                                   PS_iMeanPSFMag_qual,PS_iMeanPSFMagErr_qual,PS_iMeanKronMag_qual,PS_iMeanKronMagErr_qual,PS_iMeanApMag,PS_iMeanApMagErr,
                                   PS_zMeanPSFMag_qual,PS_zMeanPSFMagErr_qual,PS_zMeanKronMag_qual,PS_zMeanKronMagErr_qual,PS_zMeanApMag,PS_zMeanApMagErr])

    panstarrsplotdir = plotdir+"panstarrs/"
    if GENERATE:
        if pathexists(panstarrsplotdir):
            os.system("rm -rf "+panstarrsplotdir)
        if VERBOSE:
            print "cleared directory "+panstarrsplotdir
    os.system("mkdir "+panstarrsplotdir)
    if VERBOSE:
            print "creating histograms of Pan-STARRS magnitudes"

    #################################### FIND MAGNITUDE ZERO POINT ###########################################

    panstarrsplotdir_cal = plotdir+"calibration/panstarrs/"

    if pathexists(panstarrsplotdir_cal):
        os.system("rm -rf "+panstarrsplotdir_cal)
    if VERBOSE:
        print "cleared directory "+panstarrsplotdir_cal

    os.system("mkdir "+panstarrsplotdir_cal)

    pscalibrateddir = calibrateddir+"panstarrs/"
    if pathexists(pscalibrateddir):
        os.system("rm -rf "+pscalibrateddir)
    if VERBOSE:
        print "cleared directory "+pscalibrateddir
    os.system("mkdir "+pscalibrateddir)

    for folder in folders:
        panstarrsplotfolderdir   = panstarrsplotdir_cal+folder+"/"
        if pathexists(panstarrsplotfolderdir):
            os.system("rm -rf "+panstarrsplotfolderdir)
        os.system("mkdir "+panstarrsplotfolderdir)
        if VERBOSE:
            print "cleared directory "+panstarrsplotfolderdir
        #calibratedfolderdir = atlascalibrateddir+folder+"/"
        calibratedfolderdir = pscalibrateddir+folder+"/"
        if pathexists(calibratedfolderdir):
            os.system("rm -rf "+calibratedfolderdir)
            if VERBOSE:
                print "cleared directory "+calibratedfolderdir
        os.system("mkdir "+calibratedfolderdir)
        if VERBOSE:
            print "created directory "+calibratedfolderdir

    if SEXTRACTOR:
        checkdirs = [sexobjdir,sexcatdir,sexbgdir]

        # clear catalog directories
        if pathexists(sexdir):
            os.system("rm -rf "+sexdir)
            if VERBOSE:
                print "cleared directory "+sexdir
        os.system("mkdir "+sexdir)
        for dir in checkdirs:
            os.system("mkdir "+dir)
            for folder in folders:
                os.system("mkdir "+dir+folder)

        # write configuration files
        sextractor_config_name,params_name,nnw_name,conv_name,catalog_name = writesexfiles(parentdir,VERBOSE)

    if LOGOUTPUT:
        logbool=True
        logstring=[]
    else:
        logbool=False

    pool_mainzeropointcal_panstarrs = multiprocessing.Pool(ncores_panstarrs)

    if VERBOSE:
        print "\n-------------- obtaining overlap between Pan-STARRS and SExtractor ---------------"
    if LOGOUTPUT:
        logstring.append("-------------- obtaining overlap between Pan-STARRS and SExtractor ---------------")

    for folder in folders:
        panstarrsplotfolderdir   = panstarrsplotdir_cal+folder+"/"
        if VERBOSE:
            print "processing folder "+folder+"..."
        if LOGOUTPUT:
            logstring.append(" > processing folder "+str(folder)+"...")
        folderdir = findfolderdir(folder)
        if VERBOSE:
            print "using data directory "+folderdir+" for folder "+folder
        if LOGOUTPUT:
            logstring.append("using data directory "+str(folderdir)+" for folder "+str(folder))
        sciencefiles = os.listdir(folderdir)
        print sciencefiles
        if MULTIPROCESSING_CAL:
            print "Multiprocessing..."
            mainzeropointcal_mp=partial(mainzeropointcal_panstarrs,panstarrsdata_trim,folder,folderdir,logbool=logbool,CHECKZP=CHECKZP)
            logstring_zip,stackscienceimages_longexp_dict_zip,zeropoints_longexp_dict_zip,zeropoint_errs_longexp_dict_zip,airmasses_longexp_dict_zip,timeobs_longexp_dict_zip,FWHM_longexp_dict_zip,mag_sex_matches_clean_longexp_dict_zip,magerr_sex_matches_clean_longexp_dict_zip,mag_ps_matches_clean_longexp_dict_zip,magerr_ps_matches_clean_longexp_dict_zip,mag_sex_matches_longexp_dict_zip,magerr_sex_matches_longexp_dict_zip,mag_ps_matches_longexp_dict_zip,magerr_ps_matches_longexp_dict_zip,zeropoints_all_longexp_dict_zip,zeropoints_err_all_longexp_dict_zip,airmasses_all_longexp_dict_zip,timeobs_all_longexp_dict_zip,FWHM_all_longexp_dict_zip,donotstackscienceimages_longexp_dict_zip,stackscienceimages_shortexp_dict_zip,zeropoints_shortexp_dict_zip,zeropoint_errs_shortexp_dict_zip,airmasses_shortexp_dict_zip,timeobs_shortexp_dict_zip,FWHM_shortexp_dict_zip,mag_sex_matches_clean_shortexp_dict_zip,magerr_sex_matches_clean_shortexp_dict_zip,mag_ps_matches_clean_shortexp_dict_zip,magerr_ps_matches_clean_shortexp_dict_zip,mag_sex_matches_shortexp_dict_zip,magerr_sex_matches_shortexp_dict_zip,mag_ps_matches_shortexp_dict_zip,magerr_ps_matches_shortexp_dict_zip,zeropoints_all_shortexp_dict_zip,zeropoints_err_all_shortexp_dict_zip,airmasses_all_shortexp_dict_zip,timeobs_all_shortexp_dict_zip,FWHM_all_shortexp_dict_zip,donotstackscienceimages_shortexp_dict_zip=zip(*pool_mainzeropointcal_panstarrs.map(mainzeropointcal_mp,sciencefiles))

            logstring                              = collectziplists(logstring_zip)

            ######################## long exposure (120 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_longexp_dict        = collectzipdicts(stackscienceimages_longexp_dict_zip)
            zeropoints_longexp_dict                = collectzipdicts(zeropoints_longexp_dict_zip)
            zeropoint_errs_longexp_dict            = collectzipdicts(zeropoint_errs_longexp_dict_zip)
            airmasses_longexp_dict                 = collectzipdicts(airmasses_longexp_dict_zip)
            timeobs_longexp_dict                   = collectzipdicts(timeobs_longexp_dict_zip)
            FWHM_longexp_dict                      = collectzipdicts(FWHM_longexp_dict_zip)
            # matched and cleaned
            mag_sex_matches_clean_longexp_dict     = collectzipdicts(mag_sex_matches_clean_longexp_dict_zip)
            magerr_sex_matches_clean_longexp_dict  = collectzipdicts(magerr_sex_matches_clean_longexp_dict_zip)
            mag_ps_matches_clean_longexp_dict      = collectzipdicts(mag_ps_matches_clean_longexp_dict_zip)
            magerr_ps_matches_clean_longexp_dict   = collectzipdicts(magerr_ps_matches_clean_longexp_dict_zip)
            # matched
            mag_sex_matches_longexp_dict           = collectzipdicts(mag_sex_matches_longexp_dict_zip)
            magerr_sex_matches_longexp_dict        = collectzipdicts(magerr_sex_matches_longexp_dict_zip)
            mag_ps_matches_longexp_dict            = collectzipdicts(mag_ps_matches_longexp_dict_zip)
            magerr_ps_matches_longexp_dict         = collectzipdicts(magerr_ps_matches_longexp_dict_zip)
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_longexp_dict            = collectzipdicts(zeropoints_all_longexp_dict_zip)
            zeropoints_err_all_longexp_dict        = collectzipdicts(zeropoints_err_all_longexp_dict_zip)
            airmasses_all_longexp_dict             = collectzipdicts(airmasses_all_longexp_dict_zip)
            timeobs_all_longexp_dict               = collectzipdicts(timeobs_all_longexp_dict_zip)
            FWHM_all_longexp_dict                  = collectzipdicts(FWHM_all_longexp_dict_zip)
            donotstackscienceimages_longexp_dict   = collectzipdicts(donotstackscienceimages_longexp_dict_zip)

            ######################## short exposure (5 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_shortexp_dict        = collectzipdicts(stackscienceimages_shortexp_dict_zip)
            zeropoints_shortexp_dict                = collectzipdicts(zeropoints_shortexp_dict_zip)
            zeropoint_errs_shortexp_dict            = collectzipdicts(zeropoint_errs_shortexp_dict_zip)
            airmasses_shortexp_dict                 = collectzipdicts(airmasses_shortexp_dict_zip)
            timeobs_shortexp_dict                   = collectzipdicts(timeobs_shortexp_dict_zip)
            FWHM_shortexp_dict                      = collectzipdicts(FWHM_shortexp_dict_zip)
            # matched and cleaned
            mag_sex_matches_clean_shortexp_dict     = collectzipdicts(mag_sex_matches_clean_shortexp_dict_zip)
            magerr_sex_matches_clean_shortexp_dict  = collectzipdicts(magerr_sex_matches_clean_shortexp_dict_zip)
            mag_ps_matches_clean_shortexp_dict      = collectzipdicts(mag_ps_matches_clean_shortexp_dict_zip)
            magerr_ps_matches_clean_shortexp_dict   = collectzipdicts(magerr_ps_matches_clean_shortexp_dict_zip)
            # matched
            mag_sex_matches_shortexp_dict           = collectzipdicts(mag_sex_matches_shortexp_dict_zip)
            magerr_sex_matches_shortexp_dict        = collectzipdicts(magerr_sex_matches_shortexp_dict_zip)
            mag_ps_matches_shortexp_dict            = collectzipdicts(mag_ps_matches_shortexp_dict_zip)
            magerr_ps_matches_shortexp_dict         = collectzipdicts(magerr_ps_matches_shortexp_dict_zip)
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_shortexp_dict            = collectzipdicts(zeropoints_all_shortexp_dict_zip)
            zeropoints_err_all_shortexp_dict        = collectzipdicts(zeropoints_err_all_shortexp_dict_zip)
            airmasses_all_shortexp_dict             = collectzipdicts(airmasses_all_shortexp_dict_zip)
            timeobs_all_shortexp_dict               = collectzipdicts(timeobs_all_shortexp_dict_zip)
            FWHM_all_shortexp_dict                  = collectzipdicts(FWHM_all_shortexp_dict_zip)
            donotstackscienceimages_shortexp_dict   = collectzipdicts(donotstackscienceimages_shortexp_dict_zip)

        else:
            print "No multi-processing..."

            logstring                              = []

            ######################## long exposure (120 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_longexp_dict        = {}
            zeropoints_longexp_dict                = {}
            zeropoint_errs_longexp_dict            = {}
            airmasses_longexp_dict                 = {}
            timeobs_longexp_dict                   = {}
            FWHM_longexp_dict                      = {}
            # matched and cleaned
            mag_sex_matches_clean_longexp_dict     = {}
            magerr_sex_matches_clean_longexp_dict  = {}
            mag_ps_matches_clean_longexp_dict      = {}
            magerr_ps_matches_clean_longexp_dict   = {}
            # matched
            mag_sex_matches_longexp_dict           = {}
            magerr_sex_matches_longexp_dict        = {}
            mag_ps_matches_longexp_dict            = {}
            magerr_ps_matches_longexp_dict         = {}
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_longexp_dict            = {}
            zeropoints_err_all_longexp_dict        = {}
            airmasses_all_longexp_dict             = {}
            timeobs_all_longexp_dict               = {}
            FWHM_all_longexp_dict                  = {}
            donotstackscienceimages_longexp_dict   = {}

            ######################## short exposure (5 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_shortexp_dict        = {}
            zeropoints_shortexp_dict                = {}
            zeropoint_errs_shortexp_dict            = {}
            airmasses_shortexp_dict                 = {}
            timeobs_shortexp_dict                   = {}
            FWHM_shortexp_dict                      = {}
            # matched and cleaned
            mag_sex_matches_clean_shortexp_dict     = {}
            magerr_sex_matches_clean_shortexp_dict  = {}
            mag_ps_matches_clean_shortexp_dict      = {}
            magerr_ps_matches_clean_shortexp_dict   = {}
            # matched
            mag_sex_matches_shortexp_dict           = {}
            magerr_sex_matches_shortexp_dict        = {}
            mag_ps_matches_shortexp_dict            = {}
            magerr_ps_matches_shortexp_dict         = {}
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_shortexp_dict            = {}
            zeropoints_err_all_shortexp_dict        = {}
            airmasses_all_shortexp_dict             = {}
            timeobs_all_shortexp_dict               = {}
            FWHM_all_shortexp_dict                  = {}
            donotstackscienceimages_shortexp_dict   = {}

            i=0
            checkfiles = ["W3_1-S001-R001-C001-r_dupe-2-new-dsff.fts"]
            for file in sciencefiles:
                #if file in checkfiles:
                if i<1:
                    logstring_temp,stackscienceimages_longexp_dict_temp,zeropoints_longexp_dict_temp,zeropoint_errs_longexp_dict_temp,airmasses_longexp_dict_temp,timeobs_longexp_dict_temp,FWHM_longexp_dict_temp,mag_sex_matches_clean_longexp_dict_temp,magerr_sex_matches_clean_longexp_dict_temp,mag_ps_matches_clean_longexp_dict_temp,magerr_ps_matches_clean_longexp_dict_temp,mag_sex_matches_longexp_dict_temp,magerr_sex_matches_longexp_dict_temp,mag_ps_matches_longexp_dict_temp,magerr_ps_matches_longexp_dict_temp,zeropoints_all_longexp_dict_temp,zeropoints_err_all_longexp_dict_temp,airmasses_all_longexp_dict_temp,timeobs_all_longexp_dict_temp,FWHM_all_longexp_dict_temp,donotstackscienceimages_longexp_dict_temp,stackscienceimages_shortexp_dict_temp,zeropoints_shortexp_dict_temp,zeropoint_errs_shortexp_dict_temp,airmasses_shortexp_dict_temp,timeobs_shortexp_dict_temp,FWHM_shortexp_dict_temp,mag_sex_matches_clean_shortexp_dict_temp,magerr_sex_matches_clean_shortexp_dict_temp,mag_ps_matches_clean_shortexp_dict_temp,magerr_ps_matches_clean_shortexp_dict_temp,mag_sex_matches_shortexp_dict_temp,magerr_sex_matches_shortexp_dict_temp,mag_ps_matches_shortexp_dict_temp,magerr_ps_matches_shortexp_dict_temp,zeropoints_all_shortexp_dict_temp,zeropoints_err_all_shortexp_dict_temp,airmasses_all_shortexp_dict_temp,timeobs_all_shortexp_dict_temp,FWHM_all_shortexp_dict_temp,donotstackscienceimages_shortexp_dict_temp=mainzeropointcal_panstarrs(panstarrsdata_trim,folder,folderdir,file,logbool=logbool,CHECKZP=CHECKZP)

                    logstring.append(logstring_temp)

                    ######################## long exposure (120 s) ########################
                    # zero slopes and 'good' zero point values
                    stackscienceimages_longexp_dict       = append_dict2dict(stackscienceimages_longexp_dict,stackscienceimages_longexp_dict_temp)
                    zeropoints_longexp_dict               = append_dict2dict(zeropoints_longexp_dict,zeropoints_longexp_dict_temp)
                    zeropoint_errs_longexp_dict           = append_dict2dict(zeropoint_errs_longexp_dict,zeropoint_errs_longexp_dict_temp)
                    airmasses_longexp_dict                = append_dict2dict(airmasses_longexp_dict,airmasses_longexp_dict_temp)
                    timeobs_longexp_dict                  = append_dict2dict(timeobs_longexp_dict,timeobs_longexp_dict_temp)
                    FWHM_longexp_dict                     = append_dict2dict(FWHM_longexp_dict,FWHM_longexp_dict_temp)
                    # matched and cleaned
                    mag_sex_matches_clean_longexp_dict    = append_dict2dict(mag_sex_matches_clean_longexp_dict,mag_sex_matches_clean_longexp_dict_temp)
                    magerr_sex_matches_clean_longexp_dict = append_dict2dict(magerr_sex_matches_clean_longexp_dict,magerr_sex_matches_clean_longexp_dict_temp)
                    mag_ps_matches_clean_longexp_dict     = append_dict2dict(mag_ps_matches_clean_longexp_dict,mag_ps_matches_clean_longexp_dict_temp)
                    magerr_ps_matches_clean_longexp_dict  = append_dict2dict(magerr_ps_matches_clean_longexp_dict,magerr_ps_matches_clean_longexp_dict_temp)
                    # matched
                    mag_sex_matches_longexp_dict          = append_dict2dict(mag_sex_matches_longexp_dict,mag_sex_matches_longexp_dict_temp)
                    magerr_sex_matches_longexp_dict       = append_dict2dict(magerr_sex_matches_longexp_dict,magerr_sex_matches_longexp_dict_temp)
                    mag_ps_matches_longexp_dict           = append_dict2dict(mag_ps_matches_longexp_dict,mag_ps_matches_longexp_dict_temp)
                    magerr_ps_matches_longexp_dict        = append_dict2dict(magerr_ps_matches_longexp_dict,magerr_ps_matches_longexp_dict_temp)
                    # non-zero slopes and 'bad' zero-point values
                    zeropoints_all_longexp_dict           = append_dict2dict(zeropoints_all_longexp_dict,zeropoints_all_longexp_dict_temp)
                    zeropoints_err_all_longexp_dict       = append_dict2dict(zeropoints_err_all_longexp_dict,zeropoints_err_all_longexp_dict_temp)
                    airmasses_all_longexp_dict            = append_dict2dict(airmasses_all_longexp_dict,airmasses_all_longexp_dict_temp)
                    timeobs_all_longexp_dict              = append_dict2dict(timeobs_all_longexp_dict,timeobs_all_longexp_dict_temp)
                    FWHM_all_longexp_dict                 = append_dict2dict(FWHM_all_longexp_dict,FWHM_all_longexp_dict_temp)
                    donotstackscienceimages_longexp_dict  = append_dict2dict(donotstackscienceimages_longexp_dict,donotstackscienceimages_longexp_dict_temp)

                    ######################## short exposure (5 s) ########################
                    # zero slopes and 'good' zero point values
                    stackscienceimages_shortexp_dict       = append_dict2dict(stackscienceimages_shortexp_dict,stackscienceimages_shortexp_dict_temp)
                    zeropoints_shortexp_dict               = append_dict2dict(zeropoints_shortexp_dict,zeropoints_shortexp_dict_temp)
                    zeropoint_errs_shortexp_dict           = append_dict2dict(zeropoint_errs_shortexp_dict,zeropoint_errs_shortexp_dict_temp)
                    airmasses_shortexp_dict                = append_dict2dict(airmasses_shortexp_dict,airmasses_shortexp_dict_temp)
                    timeobs_shortexp_dict                  = append_dict2dict(timeobs_shortexp_dict,timeobs_shortexp_dict_temp)
                    FWHM_shortexp_dict                     = append_dict2dict(FWHM_shortexp_dict,FWHM_shortexp_dict_temp)
                    # matched and cleaned
                    mag_sex_matches_clean_shortexp_dict    = append_dict2dict(mag_sex_matches_clean_shortexp_dict,mag_sex_matches_clean_shortexp_dict_temp)
                    magerr_sex_matches_clean_shortexp_dict = append_dict2dict(magerr_sex_matches_clean_shortexp_dict,magerr_sex_matches_clean_shortexp_dict_temp)
                    mag_ps_matches_clean_shortexp_dict     = append_dict2dict(mag_ps_matches_clean_shortexp_dict,mag_ps_matches_clean_shortexp_dict_temp)
                    magerr_ps_matches_clean_shortexp_dict  = append_dict2dict(magerr_ps_matches_clean_shortexp_dict,magerr_ps_matches_clean_shortexp_dict_temp)
                    # matched
                    mag_sex_matches_shortexp_dict          = append_dict2dict(mag_sex_matches_shortexp_dict,mag_sex_matches_shortexp_dict_temp)
                    magerr_sex_matches_shortexp_dict       = append_dict2dict(magerr_sex_matches_shortexp_dict,magerr_sex_matches_shortexp_dict_temp)
                    mag_ps_matches_shortexp_dict           = append_dict2dict(mag_ps_matches_shortexp_dict,mag_ps_matches_shortexp_dict_temp)
                    magerr_ps_matches_shortexp_dict        = append_dict2dict(magerr_ps_matches_shortexp_dict,magerr_ps_matches_shortexp_dict_temp)
                    # non-zero slopes and 'bad' zero-point values
                    zeropoints_all_shortexp_dict           = append_dict2dict(zeropoints_all_shortexp_dict,zeropoints_all_shortexp_dict_temp)
                    zeropoints_err_all_shortexp_dict       = append_dict2dict(zeropoints_err_all_shortexp_dict,zeropoints_err_all_shortexp_dict_temp)
                    airmasses_all_shortexp_dict            = append_dict2dict(airmasses_all_shortexp_dict,airmasses_all_shortexp_dict_temp)
                    timeobs_all_shortexp_dict              = append_dict2dict(timeobs_all_shortexp_dict,timeobs_all_shortexp_dict_temp)
                    FWHM_all_shortexp_dict                 = append_dict2dict(FWHM_all_shortexp_dict,FWHM_all_shortexp_dict_temp)
                    donotstackscienceimages_shortexp_dict  = append_dict2dict(donotstackscienceimages_shortexp_dict,donotstackscienceimages_shortexp_dict_temp)

                    i+=1

        ######################## long exposure (120 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_longexp_dict_flat         = flattendict(stackscienceimages_longexp_dict)
        zeropoints_longexp_dict_flat                 = flattendict(zeropoints_longexp_dict)
        zeropoint_errs_longexp_dict_flat             = flattendict(zeropoint_errs_longexp_dict)
        airmasses_longexp_dict_flat                  = flattendict(airmasses_longexp_dict)
        timeobs_longexp_dict_flat                    = flattendict(timeobs_longexp_dict)
        FWHM_longexp_dict_flat                       = flattendict(FWHM_longexp_dict)
        # matched and cleaned
        mag_sex_matches_clean_longexp_dict_flat      = flattendict(mag_sex_matches_clean_longexp_dict,ARRAYS=True)
        magerr_sex_matches_clean_longexp_dict_flat   = flattendict(magerr_sex_matches_clean_longexp_dict,ARRAYS=True)
        mag_ps_matches_clean_longexp_dict_flat       = flattendict(mag_ps_matches_clean_longexp_dict,ARRAYS=True)
        magerr_ps_matches_clean_longexp_dict_flat    = flattendict(magerr_ps_matches_clean_longexp_dict,ARRAYS=True)
        # matched
        mag_sex_matches_longexp_dict_flat            = flattendict(mag_sex_matches_longexp_dict,ARRAYS=True)
        magerr_sex_matches_longexp_dict_flat         = flattendict(magerr_sex_matches_longexp_dict,ARRAYS=True)
        mag_ps_matches_longexp_dict_flat             = flattendict(mag_ps_matches_longexp_dict,ARRAYS=True)
        magerr_ps_matches_longexp_dict_flat          = flattendict(magerr_ps_matches_longexp_dict,ARRAYS=True)
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_longexp_dict_flat             = flattendict(zeropoints_all_longexp_dict)
        zeropoints_err_all_longexp_dict_flat         = flattendict(zeropoints_err_all_longexp_dict)
        airmasses_all_longexp_dict_flat              = flattendict(airmasses_all_longexp_dict)
        timeobs_all_longexp_dict_flat                = flattendict(timeobs_all_longexp_dict)
        FWHM_all_longexp_dict_flat                   = flattendict(FWHM_all_longexp_dict)
        donotstackscienceimages_longexp_dict_flat    = flattendict(donotstackscienceimages_longexp_dict)

        ######################## short exposure (5 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_shortexp_dict_flat        = flattendict(stackscienceimages_shortexp_dict)
        zeropoints_shortexp_dict_flat                = flattendict(zeropoints_shortexp_dict)
        zeropoint_errs_shortexp_dict_flat            = flattendict(zeropoint_errs_shortexp_dict)
        airmasses_shortexp_dict_flat                 = flattendict(airmasses_shortexp_dict)
        timeobs_shortexp_dict_flat                   = flattendict(timeobs_shortexp_dict)
        FWHM_shortexp_dict_flat                      = flattendict(FWHM_shortexp_dict)
        # matched and cleaned
        mag_sex_matches_clean_shortexp_dict_flat     = flattendict(mag_sex_matches_clean_shortexp_dict,ARRAYS=True)
        magerr_sex_matches_clean_shortexp_dict_flat  = flattendict(magerr_sex_matches_clean_shortexp_dict,ARRAYS=True)
        mag_ps_matches_clean_shortexp_dict_flat      = flattendict(mag_ps_matches_clean_shortexp_dict,ARRAYS=True)
        magerr_ps_matches_clean_shortexp_dict_flat   = flattendict(magerr_ps_matches_clean_shortexp_dict,ARRAYS=True)
        # matched
        mag_sex_matches_shortexp_dict_flat           = flattendict(mag_sex_matches_shortexp_dict,ARRAYS=True)
        magerr_sex_matches_shortexp_dict_flat        = flattendict(magerr_sex_matches_shortexp_dict,ARRAYS=True)
        mag_ps_matches_shortexp_dict_flat            = flattendict(mag_ps_matches_shortexp_dict,ARRAYS=True)
        magerr_ps_matches_shortexp_dict_flat         = flattendict(magerr_ps_matches_shortexp_dict,ARRAYS=True)
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_shortexp_dict_flat            = flattendict(zeropoints_all_shortexp_dict)
        zeropoints_err_all_shortexp_dict_flat        = flattendict(zeropoints_err_all_shortexp_dict)
        airmasses_all_shortexp_dict_flat             = flattendict(airmasses_all_shortexp_dict)
        timeobs_all_shortexp_dict_flat               = flattendict(timeobs_all_shortexp_dict)
        FWHM_all_shortexp_dict_flat                  = flattendict(FWHM_all_shortexp_dict)
        donotstackscienceimages_shortexp_dict_flat   = flattendict(donotstackscienceimages_shortexp_dict)

        # store dictionary of files to stack
        fwritedict_pickle(stackscienceimages_longexp_dict_flat,"stackscienceimages_longexp_dict_pickle.txt")
        fwritedict_pickle(stackscienceimages_shortexp_dict_flat,"stackscienceimages_shortexp_dict_pickle.txt")
        fwritedict_text(stackscienceimages_longexp_dict_flat,"stackscienceimages_longexp_dict.txt")
        fwritedict_text(stackscienceimages_shortexp_dict_flat,"stackscienceimages_shortexp_dict.txt")

        if VERBOSE:
            print "Plotting all zero-point data for each passband..."

        # plot zeropoint data for each passband
        plotzp_all(mag_sex_matches_clean_longexp_dict_flat,magerr_sex_matches_clean_longexp_dict_flat,mag_ps_matches_clean_longexp_dict_flat,magerr_ps_matches_clean_longexp_dict_flat,mag_sex_matches_longexp_dict_flat,magerr_sex_matches_longexp_dict_flat,mag_ps_matches_longexp_dict_flat,magerr_ps_matches_longexp_dict_flat,panstarrsplotfolderdir,"LONG",MULTIPROCESSING_CAL,CHECKZP)
        plotzp_all(mag_sex_matches_clean_shortexp_dict_flat,magerr_sex_matches_clean_shortexp_dict_flat,mag_ps_matches_clean_shortexp_dict_flat,magerr_ps_matches_clean_shortexp_dict_flat,mag_sex_matches_shortexp_dict_flat,magerr_sex_matches_shortexp_dict_flat,mag_ps_matches_shortexp_dict_flat,magerr_ps_matches_shortexp_dict_flat,panstarrsplotfolderdir,"SHORT",MULTIPROCESSING_CAL,CHECKZP)

        if VERBOSE:
            print "Plotting zero-point correlations for each passband..."

        # plot zeropoint correlations
        plotzpcorrs(zeropoints_longexp_dict_flat,zeropoint_errs_longexp_dict_flat,zeropoints_all_longexp_dict_flat,zeropoints_err_all_longexp_dict_flat,airmasses_longexp_dict_flat,airmasses_all_longexp_dict_flat,timeobs_longexp_dict_flat,timeobs_all_longexp_dict_flat,FWHM_longexp_dict_flat,FWHM_all_longexp_dict_flat,"LONG",panstarrsplotfolderdir,CHECKZP)
        plotzpcorrs(zeropoints_shortexp_dict_flat,zeropoint_errs_shortexp_dict_flat,zeropoints_all_shortexp_dict_flat,zeropoints_err_all_shortexp_dict_flat,airmasses_shortexp_dict_flat,airmasses_all_shortexp_dict_flat,timeobs_shortexp_dict_flat,timeobs_all_shortexp_dict_flat,FWHM_shortexp_dict_flat,FWHM_all_shortexp_dict_flat,"SHORT",panstarrsplotfolderdir,CHECKZP)

        # clear variables from memory (running the initial calibration and the second calibration step back-to-back will double their memory and raise a memory error)
        if VERBOSE:
            print "clearing calibration vairables to minimize memory..."

        ######################## long exposure (120 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_longexp_dict = {}             # non-flat dictionaries
        zeropoints_longexp_dict = {}                     # non-flat dictionaries
        zeropoint_errs_longexp_dict = {}                 # non-flat dictionaries
        airmasses_longexp_dict_flat = {}                 # for ZP correlation plots
        airmasses_longexp_dict = {}                      # for ZP correlation plots
        timeobs_longexp_dict_flat = {}                   # for ZP correlation plots
        timeobs_longexp_dict = {}                        # for ZP correlation plots
        FWHM_longexp_dict_flat = {}                      # for ZP correlation plots
        FWHM_longexp_dict = {}                           # for ZP correlation plots
        # matched and cleaned
        mag_sex_matches_clean_longexp_dict_flat = {}     # for ZP plots
        mag_sex_matches_clean_longexp_dict = {}          # for ZP plots
        magerr_sex_matches_clean_longexp_dict_flat = {}  # for ZP plots
        magerr_sex_matches_clean_longexp_dict = {}       # for ZP plots
        mag_ps_matches_clean_longexp_dict_flat = {}      # for ZP plots
        mag_ps_matches_clean_longexp_dict = {}           # for ZP plots
        magerr_ps_matches_clean_longexp_dict_flat = {}   # for ZP plots
        magerr_ps_matches_clean_longexp_dict = {}        # for ZP plots
        # matched
        mag_sex_matches_longexp_dict_flat = {}           # for ZP plots
        mag_sex_matches_longexp_dict = {}                # for ZP plots
        magerr_sex_matches_longexp_dict_flat = {}        # for ZP plots
        magerr_sex_matches_longexp_dict = {}             # for ZP plots
        mag_ps_matches_longexp_dict_flat = {}            # for ZP plots
        mag_ps_matches_longexp_dict = {}                 # for ZP plots
        magerr_ps_matches_longexp_dict_flat = {}         # for ZP plots
        magerr_ps_matches_longexp_dict = {}              # for ZP plots
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_longexp_dict_flat = {}            # for ZP plots
        zeropoints_all_longexp_dict = {}                 # for ZP plots
        zeropoints_err_all_longexp_dict_flat = {}        # for ZP plots
        zeropoints_err_all_longexp_dict = {}             # for ZP plots
        airmasses_all_longexp_dict_flat = {}             # for ZP correlation plots
        airmasses_all_longexp_dict = {}                  # for ZP correlation plots
        timeobs_all_longexp_dict_flat = {}               # for ZP correlation plots
        timeobs_all_longexp_dict = {}                    # for ZP correlation plots
        FWHM_all_longexp_dict_flat = {}                  # for ZP correlation plots
        FWHM_all_longexp_dict = {}                       # for ZP correlation plots
        donotstackscienceimages_longexp_dict = {}        # non-flat dictionary
        ######################## short exposure (5 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_shortexp_dict = {}            # non-flat dictionaries
        zeropoints_shortexp_dict = {}                    # non-flat dictionaries
        zeropoint_errs_shortexp_dict = {}                # non-flat dictionaries
        airmasses_shortexp_dict_flat = {}                # for ZP correlation plots
        airmasses_shortexp_dict = {}                     # for ZP correlation plots
        timeobs_shortexp_dict_flat = {}                  # for ZP correlation plots
        timeobs_shortexp_dict = {}                       # for ZP correlation plots
        FWHM_shortexp_dict_flat = {}                     # for ZP correlation plots
        FWHM_shortexp_dict = {}                          # for ZP correlation plots
        # matched and cleaned
        mag_sex_matches_clean_shortexp_dict_flat = {}    # for ZP plots
        mag_sex_matches_clean_shortexp_dict = {}         # for ZP plots
        magerr_sex_matches_clean_shortexp_dict_flat = {} # for ZP plots
        magerr_sex_matches_clean_shortexp_dict = {}      # for ZP plots
        mag_ps_matches_clean_shortexp_dict_flat = {}     # for ZP plots
        mag_ps_matches_clean_shortexp_dict = {}          # for ZP plots
        magerr_ps_matches_clean_shortexp_dict_flat = {}  # for ZP plots
        magerr_ps_matches_clean_shortexp_dict = {}       # for ZP plots
        # matched
        mag_sex_matches_shortexp_dict_flat = {}          # for ZP plots
        mag_sex_matches_shortexp_dict = {}               # for ZP plots
        magerr_sex_matches_shortexp_dict_flat = {}       # for ZP plots
        magerr_sex_matches_shortexp_dict = {}            # for ZP plots
        mag_ps_matches_shortexp_dict_flat = {}           # for ZP plots
        mag_ps_matches_shortexp_dict = {}                # for ZP plots
        magerr_ps_matches_shortexp_dict_flat = {}        # for ZP plots
        magerr_ps_matches_shortexp_dict = {}             # for ZP plots
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_shortexp_dict_flat = {}           # for ZP plots
        zeropoints_all_shortexp_dict = {}                # for ZP plots
        zeropoints_err_all_shortexp_dict_flat = {}       # for ZP plots
        zeropoints_err_all_shortexp_dict = {}            # for ZP plots
        airmasses_all_shortexp_dict_flat = {}            # for ZP correlation plots
        airmasses_all_shortexp_dict = {}                 # for ZP correlation plots
        timeobs_all_shortexp_dict_flat = {}              # for ZP correlation plots
        timeobs_all_shortexp_dict = {}                   # for ZP correlation plots
        FWHM_all_shortexp_dict_flat = {}                 # for ZP correlation plots
        FWHM_all_shortexp_dict = {}                      # for ZP correlation plots
        donotstackscienceimages_shortexp_dict = {}       # non-flat dictionary

        print "\nscience images to be stacked (long exp):"
        print stackscienceimages_longexp_dict_flat
        print "\nbad science images (long exposure):"
        print donotstackscienceimages_longexp_dict_flat
        print "\nscience images to be stacked (short exp):"
        print stackscienceimages_shortexp_dict_flat
        print "\nbad science images (short exposure):"
        print donotstackscienceimages_shortexp_dict_flat
        #print "\nlogstring:"
        #print logstring

        if LOGOUTPUT:
            # write calibration log
            if VERBOSE:
                print "writing zeropoint log"
            f=open(logdir+"zeropoint.log","w+")
            for item in logstring:
                f.write("%s\n" % item)
            f.close()

else:
    if VERBOSE:
        print "\n------------------------- skipping Pan-STARRS calilbration -----------------------"

#########################################################################################################
#                                                                                                       #
#                                        ATLAS CALIBRATION                                              #
#                                                                                                       #
#########################################################################################################

#########################################################################################################
#                                                                                                       #
#                                     ZEROPOINT CALIBRATION                                             #
#                                                                                                       #
#########################################################################################################

######################################## ZEROPOINT FUNCTION #############################################

atlaspassbands = ["g","r","i","z"]

def mainzeropointcal_atlas(atlasdata,folder,folderdir,file,logbool=False,CHECKZP=False,FINALCHECKZP=False):
    #
    # Measures the magnitude zero point by comparing Source Extractor measurements to ATLAS. Uses the slope of
    # (ATLAS mag - SExtractor mag) vs. ATLAS mag to determine if the image will be included; a non-zero slope
    # within 3-sigma uncertainties will be returned in the 'badcalfiles' dictionary and not used any further.
    #

    logstring                                = []

    ######################## long exposure (120 s) ########################
    # zero slopes and 'good' zero point values
    stackscienceimages_longexp_dict          = {}
    zeropoints_longexp_dict                  = {}
    zeropoint_errs_longexp_dict              = {}
    airmasses_longexp_dict                   = {}
    timeobs_longexp_dict                     = {}
    FWHM_longexp_dict                        = {}
    # matched and cleaned
    mag_sex_matches_clean_longexp_dict       = {}
    magerr_sex_matches_clean_longexp_dict    = {}
    mag_atlas_matches_clean_longexp_dict     = {}
    magerr_atlas_matches_clean_longexp_dict  = {}
    # matched
    mag_sex_matches_longexp_dict             = {}
    magerr_sex_matches_longexp_dict          = {}
    mag_atlas_matches_longexp_dict           = {}
    magerr_atlas_matches_longexp_dict        = {}
    # non-zero slopes and 'bad' zero-point values
    zeropoints_all_longexp_dict              = {}
    zeropoints_err_all_longexp_dict          = {}
    airmasses_all_longexp_dict               = {}
    timeobs_all_longexp_dict                 = {}
    FWHM_all_longexp_dict                    = {}
    donotstackscienceimages_longexp_dict     = {}

    ######################## short exposure (5 s) ########################
    # zero slopes and 'good' zero point values
    stackscienceimages_shortexp_dict         = {}
    zeropoints_shortexp_dict                 = {}
    zeropoint_errs_shortexp_dict             = {}
    airmasses_shortexp_dict                  = {}
    timeobs_shortexp_dict                    = {}
    FWHM_shortexp_dict                       = {}
    # matched and cleaned
    mag_sex_matches_clean_shortexp_dict      = {}
    magerr_sex_matches_clean_shortexp_dict   = {}
    mag_atlas_matches_clean_shortexp_dict    = {}
    magerr_atlas_matches_clean_shortexp_dict = {}
    # matched
    mag_sex_matches_shortexp_dict            = {}
    magerr_sex_matches_shortexp_dict         = {}
    mag_atlas_matches_shortexp_dict          = {}
    magerr_atlas_matches_shortexp_dict       = {}
    # non-zero slopes and 'bad' zero-point values
    zeropoints_all_shortexp_dict             = {}
    zeropoints_err_all_shortexp_dict         = {}
    airmasses_all_shortexp_dict              = {}
    timeobs_all_shortexp_dict                = {}
    FWHM_all_shortexp_dict                   = {}
    donotstackscienceimages_shortexp_dict    = {}

    sciencedata   = fits.getdata(folderdir+file)
    scienceheader = fits.getheader(folderdir+file)
    passband      = scienceheader["FILTER"]
    airmass       = scienceheader["AIRMASS"]
    exptime       = scienceheader["EXPTIME"]
    timeobs       = scienceheader["TIME-OBS"]
    timeobs_24h   = ftimeobs_to24h(timeobs)
    FWHM          = scienceheader["FWHM"]
    object        = re.split("-|\.",file)[0]
    if passband in atlaspassbands:
        if VERBOSE:
            print "    > processing file: "+file
        if LOGOUTPUT:
            logstring.append("    > processing file: "+str(file))

        # run Source Extractor if indicated otherwise open already existing file
        if (CHECKZP==False) and (FINALCHECKZP==False):
            # initial calibration
            sexcat = sexcall(file,folder,folderdir,sexobjdir,sexcatdir,sexbgdir,sexbool=SEXTRACTOR,CHECKZP=CHECKZP,FINALCHECKZP=FINALCHECKZP)
        elif CHECKZP==True:
            # first calibration check
            sexcat = sexcall(file,folder,folderdir,sexobjcheckzpdir,sexcatcheckzpdir,sexbgcheckzpdir,sexbool=SEXTRACTOR,CHECKZP=CHECKZP,FINALCHECKZP=FINALCHECKZP)
        elif FINALCHECKZP==True:
            # final calibration check
            sexcat = sexcall(file,folder,folderdir,sexobjfinalcheckzpdir,sexcatfinalcheckzpdir,sexbgfinalcheckzpdir,sexbool=SEXTRACTOR,CHECKZP=CHECKZP,FINALCHECKZP=FINALCHECKZP)

        # find ATLAS overlap with SExtractor
        ZEROSLOPE,zeropoint,zeropoint_error,GOODZP,CHECKOBSTIME,mag_sex_matches_clean,magerr_sex_matches_clean,mag_atlas_matches_clean,magerr_atlas_matches_clean,mag_sex_matches,magerr_sex_matches,mag_atlas_matches,magerr_atlas_matches = magcalibration(file,folder,folderdir,atlasdata,sexcat,VERBOSE,GENERATE,CHECKZP,FINALCHECKZP)

        ######################## long exposure (120 s) ########################
        if exptime==120.:
            # collect files with zero slopes and 'good' zero point values
            if (ZEROSLOPE==True) and (GOODZP==True) and (CHECKOBSTIME==True):
                # zero slopes and 'good' zero point values
                stackscienceimages_longexp_dict          = appenditemtodict(folderdir+file,object,stackscienceimages_longexp_dict)
                zeropoints_longexp_dict                  = appenditemtodict(zeropoint,passband,zeropoints_longexp_dict)
                zeropoint_errs_longexp_dict              = appenditemtodict(zeropoint_error,passband,zeropoint_errs_longexp_dict)
                airmasses_longexp_dict                   = appenditemtodict(airmass,passband,airmasses_longexp_dict)
                timeobs_longexp_dict                     = appenditemtodict(timeobs_24h,passband,timeobs_longexp_dict)
                FWHM_longexp_dict                        = appenditemtodict(FWHM,passband,FWHM_longexp_dict)
                # matched and cleaned
                mag_sex_matches_clean_longexp_dict       = appenditemtodict(mag_sex_matches_clean,passband,mag_sex_matches_clean_longexp_dict)
                magerr_sex_matches_clean_longexp_dict    = appenditemtodict(magerr_sex_matches_clean,passband,magerr_sex_matches_clean_longexp_dict)
                mag_atlas_matches_clean_longexp_dict     = appenditemtodict(mag_atlas_matches_clean,passband,mag_atlas_matches_clean_longexp_dict)
                magerr_atlas_matches_clean_longexp_dict  = appenditemtodict(magerr_atlas_matches_clean,passband,magerr_atlas_matches_clean_longexp_dict)
                # matched
                mag_sex_matches_longexp_dict             = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_longexp_dict)
                magerr_sex_matches_longexp_dict          = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_longexp_dict)
                mag_atlas_matches_longexp_dict           = appenditemtodict(mag_atlas_matches,passband,mag_atlas_matches_longexp_dict)
                magerr_atlas_matches_longexp_dict        = appenditemtodict(magerr_atlas_matches,passband,magerr_atlas_matches_longexp_dict)
                if (CHECKZP==False) and (FINALCHECKZP==False):
                    # add zero-point to FITS header
                    fits.setval(folderdir+file, "ZP", value=zeropoint)
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=True)
                    # write calibrated image
                    calibratedfolderdir    = atlascalibrateddir+folder+"/"
                    newfilename            = file.split(".fts")[0]+"_calibrated.fts"
                    #sciencedata_cal       = fZPcorrADUs(sciencedata,zeropoint)
                    if VERBOSE:
                        print "writing calibrated image: "+calibratedfolderdir+newfilename
                    fits.writeto(calibratedfolderdir+newfilename,sciencedata,scienceheader,overwrite=True)
                elif FINALCHECKZP==True:
                    fits.setval(folderdir+file, "ZP", value=zeropoint)

            # collect files with non-zero slopes and 'bad' zero-point values
            else:
                donotstackscienceimages_longexp_dict     = appenditemtodict(folderdir+file,object,donotstackscienceimages_longexp_dict)
                if zeropoint!=None:
                    zeropoints_all_longexp_dict          = appenditemtodict(zeropoint,passband,zeropoints_all_longexp_dict)
                    zeropoints_err_all_longexp_dict      = appenditemtodict(zeropoint_error,passband,zeropoints_err_all_longexp_dict)
                    airmasses_all_longexp_dict           = appenditemtodict(airmass,passband,airmasses_all_longexp_dict)
                    timeobs_all_longexp_dict             = appenditemtodict(timeobs_24h,passband,timeobs_all_longexp_dict)
                    FWHM_all_longexp_dict                = appenditemtodict(FWHM,passband,FWHM_all_longexp_dict)
                    # matched
                    mag_sex_matches_longexp_dict         = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_longexp_dict)
                    magerr_sex_matches_longexp_dict      = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_longexp_dict)
                    mag_atlas_matches_longexp_dict       = appenditemtodict(mag_atlas_matches,passband,mag_atlas_matches_longexp_dict)
                    magerr_atlas_matches_longexp_dict    = appenditemtodict(magerr_atlas_matches,passband,magerr_atlas_matches_longexp_dict)
                if (CHECKZP==False) and (FINALCHECKZP==False):
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=False)

        ######################## short exposure (5 s) ########################
        elif exptime==5.:
            # collect files with zero slopes and 'good' zero point values
            if (ZEROSLOPE==True) and (GOODZP==True) and (CHECKOBSTIME==True):
                # zero slopes and 'good' zero point values
                stackscienceimages_shortexp_dict         = appenditemtodict(folderdir+file,object,stackscienceimages_shortexp_dict)
                zeropoints_shortexp_dict                 = appenditemtodict(zeropoint,passband,zeropoints_shortexp_dict)
                zeropoint_errs_shortexp_dict             = appenditemtodict(zeropoint_error,passband,zeropoint_errs_shortexp_dict)
                airmasses_shortexp_dict                  = appenditemtodict(airmass,passband,airmasses_shortexp_dict)
                timeobs_shortexp_dict                    = appenditemtodict(timeobs_24h,passband,timeobs_shortexp_dict)
                FWHM_shortexp_dict                       = appenditemtodict(FWHM,passband,FWHM_shortexp_dict)
                # matched and cleaned
                mag_sex_matches_clean_shortexp_dict      = appenditemtodict(mag_sex_matches_clean,passband,mag_sex_matches_clean_shortexp_dict)
                magerr_sex_matches_clean_shortexp_dict   = appenditemtodict(magerr_sex_matches_clean,passband,magerr_sex_matches_clean_shortexp_dict)
                mag_atlas_matches_clean_shortexp_dict    = appenditemtodict(mag_atlas_matches_clean,passband,mag_atlas_matches_clean_shortexp_dict)
                magerr_atlas_matches_clean_shortexp_dict = appenditemtodict(magerr_atlas_matches_clean,passband,magerr_atlas_matches_clean_shortexp_dict)
                # matched
                mag_sex_matches_shortexp_dict            = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_shortexp_dict)
                magerr_sex_matches_shortexp_dict         = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_shortexp_dict)
                mag_atlas_matches_shortexp_dict          = appenditemtodict(mag_atlas_matches,passband,mag_atlas_matches_shortexp_dict)
                magerr_atlas_matches_shortexp_dict       = appenditemtodict(magerr_atlas_matches,passband,magerr_atlas_matches_shortexp_dict)
                if (CHECKZP==False) and (FINALCHECKZP==False):
                    # add zero-point to FITS header
                    fits.setval(folderdir+file, "ZP", value=zeropoint)
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=True)
                    # write calibrated image
                    calibratedfolderdir    = atlascalibrateddir+folder+"/"
                    newfilename            = file.split(".fts")[0]+"_calibrated.fts"
                    #sciencedata_cal       = fZPcorrADUs(sciencedata,zeropoint)
                    if VERBOSE:
                        print "writing calibrated image: "+calibratedfolderdir+newfilename
                    fits.writeto(calibratedfolderdir+newfilename,sciencedata,scienceheader,overwrite=True)
                elif FINALCHECKZP==True:
                    fits.setval(folderdir+file, "ZP", value=zeropoint)
            else:
                # collect files with non-zero slopes and 'bad' zero-point values
                donotstackscienceimages_shortexp_dict    = appenditemtodict(folderdir+file,object,donotstackscienceimages_shortexp_dict)
                if zeropoint!=None:
                    zeropoints_all_shortexp_dict         = appenditemtodict(zeropoint,passband,zeropoints_all_shortexp_dict)
                    zeropoints_err_all_shortexp_dict     = appenditemtodict(zeropoint_error,passband,zeropoints_err_all_shortexp_dict)
                    airmasses_all_shortexp_dict          = appenditemtodict(airmass,passband,airmasses_all_shortexp_dict)
                    timeobs_all_shortexp_dict            = appenditemtodict(timeobs_24h,passband,timeobs_all_shortexp_dict)
                    FWHM_all_shortexp_dict               = appenditemtodict(FWHM,passband,FWHM_all_shortexp_dict)
                    # matched
                    mag_sex_matches_shortexp_dict        = appenditemtodict(mag_sex_matches,passband,mag_sex_matches_shortexp_dict)
                    magerr_sex_matches_shortexp_dict     = appenditemtodict(magerr_sex_matches,passband,magerr_sex_matches_shortexp_dict)
                    mag_atlas_matches_shortexp_dict      = appenditemtodict(mag_atlas_matches,passband,mag_atlas_matches_shortexp_dict)
                    magerr_atlas_matches_shortexp_dict   = appenditemtodict(magerr_atlas_matches,passband,magerr_atlas_matches_shortexp_dict)
                if (CHECKZP==False) and (FINALCHECKZP==False):
                    # add stacking boolean
                    fits.setval(folderdir+file, "stack", value=False)

    return logstring,stackscienceimages_longexp_dict,zeropoints_longexp_dict,zeropoint_errs_longexp_dict,airmasses_longexp_dict,timeobs_longexp_dict,FWHM_longexp_dict,mag_sex_matches_clean_longexp_dict,magerr_sex_matches_clean_longexp_dict,mag_atlas_matches_clean_longexp_dict,magerr_atlas_matches_clean_longexp_dict,mag_sex_matches_longexp_dict,magerr_sex_matches_longexp_dict,mag_atlas_matches_longexp_dict,magerr_atlas_matches_longexp_dict,zeropoints_all_longexp_dict,zeropoints_err_all_longexp_dict,airmasses_all_longexp_dict,timeobs_all_longexp_dict,FWHM_all_longexp_dict,donotstackscienceimages_longexp_dict,stackscienceimages_shortexp_dict,zeropoints_shortexp_dict,zeropoint_errs_shortexp_dict,airmasses_shortexp_dict,timeobs_shortexp_dict,FWHM_shortexp_dict,mag_sex_matches_clean_shortexp_dict,magerr_sex_matches_clean_shortexp_dict,mag_atlas_matches_clean_shortexp_dict,magerr_atlas_matches_clean_shortexp_dict,mag_sex_matches_shortexp_dict,magerr_sex_matches_shortexp_dict,mag_atlas_matches_shortexp_dict,magerr_atlas_matches_shortexp_dict,zeropoints_all_shortexp_dict,zeropoints_err_all_shortexp_dict,airmasses_all_shortexp_dict,timeobs_all_shortexp_dict,FWHM_all_shortexp_dict,donotstackscienceimages_shortexp_dict

if CALIBRATE and ATLAS:
    CHECKZP=False

    # read in catalogue data
    atlascatf = atlasdir+atlascat

    if VERBOSE:
        print "\n------------------------ proceeding with ATLAS calilbration ----------------------"
    print "Reading in ATLAS catalogue..."
    hdul_atlas       = fits.open(atlascatf)
    atlasdata        = hdul_atlas[1].data

    ATLAS_objid            = atlasdata["objid"]    # Object ID                                                                                           [none]
    ATLAS_RA               = atlasdata["RA"]       # Right ascension from Gaia DR2, J2000, epoch 2015.5                                                  [deg]
    ATLAS_DEC              = atlasdata["DEC"]      # Declination from Gaia DR2, J2000, epoch 2015.5                                                      [deg]
    ATLAS_plx              = atlasdata["plx"]      # Parallax from Gaia DR2                                                                              [mas]
    ATLAS_dplx             = atlasdata["dplx"]     # Parallax uncertainty from Gaia DR2                                                                  [mas]
    ATLAS_pmra             = atlasdata["pmra"]     # Proper motion in right ascension from Gaia DR2                                                      [mas/yr]
    ATLAS_dpmra            = atlasdata["dpmra"]    # Proper motion uncertainty in right ascension                                                        [mas/yr]
    ATLAS_pmdec            = atlasdata["pmdec"]    # Proper motion in declination from Gaia DR2                                                          [mas/yr]
    ATLAS_dpmdec           = atlasdata["dpmdec"]   # Proper motion uncertainty in declination                                                            [mas/yr]
    ATLAS_Gaia             = atlasdata["Gaia"]     # Gaia G magnitude                                                                                    [mag]
    ATLAS_dGaia            = atlasdata["dGaia"]    # Gaia G magnitude uncertainty                                                                        [mag]
    ATLAS_BP               = atlasdata["BP"]       # Gaia G_bp magnitude                                                                                 [mag]
    ATLAS_dBP              = atlasdata["dBP"]      # Gaia G_bp magnitude uncertainty                                                                     [mag]
    ATLAS_RP               = atlasdata["RP"]       # Gaia G_rp magnitude                                                                                 [mag]
    ATLAS_dRP              = atlasdata["dRP"]      # Gaia G_rp magnitude uncertainty                                                                     [mag]
    ATLAS_Teff             = atlasdata["Teff"]     # Gaia stellar effective temperature                                                                  [K]
    ATLAS_AGaia            = atlasdata["AGaia"]    # Gaia estimate of G-band extinction for this star                                                    [mag]
    ATLAS_dupvar           = atlasdata["dupvar"]   # Gaia variability and duplicate flags, 0/1/2 for "CONSTANT"/"VARIABLE"/"NOT AVAILABLE" + 4*DUPLICATE [none]
    ATLAS_Ag               = atlasdata["Ag"]       # SFD estimate of total g-band extinction                                                             [mag]
    ATLAS_rp1              = atlasdata["rp1"]      # Radius where cummulative G flux exceeds 0.1 x this star                                             [arcsec]
    ATLAS_r1               = atlasdata["r1"]       # Radius where cummulative G flux exceeds 1.0 x this star                                             [arcsec]
    ATLAS_r10              = atlasdata["r10"]      # Radius where cummulative G flux exceeds 10.0 x this star                                            [arcsec]
    ATLAS_g                = atlasdata["g"]        # PanSTARRS g magnitude                                                                               [mag]
    ATLAS_dg               = atlasdata["dg"]       # PanSTARRS g magnitude uncertainty                                                                   [mag]
    ATLAS_gchi             = atlasdata["gchi"]     # chi^2 / DOF for contributors                                                                        [none]
    ATLAS_gcontrib         = atlasdata["gcontrib"] # Bitmap of conributing catalogs to g                                                                 [none]
    ATLAS_r                = atlasdata["r"]        # PanSTARRS r magnitude                                                                               [mag]
    ATLAS_dr               = atlasdata["dr"]       # PanSTARRS r magnitude uncertainty                                                                   [mag]
    ATLAS_rchi             = atlasdata["rchi"]     # chi^2 / DOF for contributors                                                                        [none]
    ATLAS_rcontrib         = atlasdata["rcontrib"] # Bitmap of conributing catalogs to r                                                                 [none]
    ATLAS_i                = atlasdata["i"]        # PanSTARRS i magnitude                                                                               [mag]
    ATLAS_di               = atlasdata["di"]       # PanSTARRS i magnitude uncertainty                                                                   [mag]
    ATLAS_ichi             = atlasdata["ichi"]     # chi^2 / DOF for contributors                                                                        [none]
    ATLAS_icontrib         = atlasdata["icontrib"] # Bitmap of conributing catalogs to i                                                                 [none]
    ATLAS_z                = atlasdata["z"]        # PanSTARRS z magnitude                                                                               [mag]
    ATLAS_dz               = atlasdata["dz"]       # PanSTARRS z magnitude uncertainty                                                                   [mag]
    ATLAS_zchi             = atlasdata["zchi"]     # chi^2 / DOF for contributors                                                                        [none]
    ATLAS_zcontrib         = atlasdata["zcontrib"] # Bitmap of conributing catalogs to z                                                                 [none]
    ATLAS_nstat            = atlasdata["nstat"]    # Count of griz outliers rejected                                                                     [none]
    ATLAS_J                = atlasdata["J"]        # 2MASS J magnitude                                                                                   [mag]
    ATLAS_dJ               = atlasdata["dJ"]       # 2MASS J magnitude uncertainty                                                                       [mag]
    ATLAS_H                = atlasdata["H"]        # 2MASS H magnitude                                                                                   [mag]
    ATLAS_dH               = atlasdata["dH"]       # 2MASS H magnitude uncertainty                                                                       [mag]
    ATLAS_K                = atlasdata["K"]        # 2MASS K magnitude                                                                                   [mag]
    ATLAS_dK               = atlasdata["dK"]       # 2MASS K magnitude uncertainty                                                                       [mag]

    # trimmed version of ATLAS data
    atlastrimcatf = atlasdir+atlastrimcat
    t = Table([ATLAS_RA,ATLAS_DEC,ATLAS_g,ATLAS_dg,ATLAS_r,ATLAS_dr,ATLAS_i,ATLAS_di,ATLAS_z,ATLAS_dz],names=("RA","DEC","g","dg","r","dr","i","di","z","dz"))
    t.write(atlastrimcatf,format="fits",overwrite=True)

    hdul_atlastrim       = fits.open(atlastrimcatf)
    atlasdata_trim       = hdul_atlastrim[1].data

    #################################### FIND MAGNITUDE ZERO POINT ###########################################

    atlasplotdir_cal  = plotdir+"calibration/atlas/"
    if pathexists(atlasplotdir_cal):
        os.system("rm -rf "+atlasplotdir_cal)
    if VERBOSE:
        print "cleared directory "+atlasplotdir_cal
    os.system("mkdir "+atlasplotdir_cal)

    atlascalibrateddir = calibrateddir+"atlas/"
    if pathexists(atlascalibrateddir):
        os.system("rm -rf "+atlascalibrateddir)
    if VERBOSE:
        print "cleared directory "+atlascalibrateddir
    os.system("mkdir "+atlascalibrateddir)

    for folder in folders:
        atlasplotfolderdir = atlasplotdir_cal+folder+"/"
        if pathexists(atlasplotfolderdir):
            os.system("rm -rf "+atlasplotfolderdir)
        os.system("mkdir "+atlasplotfolderdir)
        if VERBOSE:
            print "cleared directory "+atlasplotfolderdir
        calibratedfolderdir = atlascalibrateddir+folder+"/"
        if pathexists(calibratedfolderdir):
            os.system("rm -rf "+calibratedfolderdir)
            if VERBOSE:
                print "cleared directory "+calibratedfolderdir
        os.system("mkdir "+calibratedfolderdir)
        if VERBOSE:
            print "created directory "+calibratedfolderdir

    if SEXTRACTOR:
        checkdirs = [sexobjdir,sexcatdir,sexbgdir]

        # clear catalog directories
        if pathexists(sexdir):
            os.system("rm -rf "+sexdir)
            if VERBOSE:
                print "cleared directory "+sexdir
        os.system("mkdir "+sexdir)
        for dir in checkdirs:
            os.system("mkdir "+dir)
            for folder in folders:
                os.system("mkdir "+dir+folder)

        # write configuration files
        #sextractor_config_name,params_name,nnw_name,conv_name,catalog_name = writesexfiles(parentdir,VERBOSE)

    if LOGOUTPUT:
        logbool=True
        logstring=[]
    else:
        logbool=False

    pool_mainzeropointcal_atlas = multiprocessing.Pool(ncores_atlas)

    if VERBOSE:
        print "\n----------------- obtaining overlap between ATLAS and SExtractor -----------------"
    if LOGOUTPUT:
        logstring.append("----------------- obtaining overlap between ATLAS and SExtractor -----------------")

    # collect for all folders; flatten later
    stackscienceimages_all_longexp_dict_temp  = {}
    stackscienceimages_all_shortexp_dict_temp = {}
    for folder in folders:
        atlasplotfolderdir    = atlasplotdir_cal+folder+"/"
        if VERBOSE:
            print "processing folder "+folder+"..."
        if LOGOUTPUT:
            logstring.append(" > processing folder "+str(folder)+"...")
        folderdir = findfolderdir(folder)
        if VERBOSE:
            print "using data directory "+folderdir+" for folder "+folder
        if LOGOUTPUT:
            logstring.append("using data directory "+str(folderdir)+" for folder "+str(folder))
        sciencefiles = os.listdir(folderdir)
        print sciencefiles
        if MULTIPROCESSING_CAL:
            print "Multiprocessing..."
            mainzeropointcal_mp=partial(mainzeropointcal_atlas,atlasdata_trim,folder,folderdir,logbool=logbool,CHECKZP=CHECKZP)
            logstring_zip,stackscienceimages_longexp_dict_zip,zeropoints_longexp_dict_zip,zeropoint_errs_longexp_dict_zip,airmasses_longexp_dict_zip,timeobs_longexp_dict_zip,FWHM_longexp_dict_zip,mag_sex_matches_clean_longexp_dict_zip,magerr_sex_matches_clean_longexp_dict_zip,mag_atlas_matches_clean_longexp_dict_zip,magerr_atlas_matches_clean_longexp_dict_zip,mag_sex_matches_longexp_dict_zip,magerr_sex_matches_longexp_dict_zip,mag_atlas_matches_longexp_dict_zip,magerr_atlas_matches_longexp_dict_zip,zeropoints_all_longexp_dict_zip,zeropoints_err_all_longexp_dict_zip,airmasses_all_longexp_dict_zip,timeobs_all_longexp_dict_zip,FWHM_all_longexp_dict_zip,donotstackscienceimages_longexp_dict_zip,stackscienceimages_shortexp_dict_zip,zeropoints_shortexp_dict_zip,zeropoint_errs_shortexp_dict_zip,airmasses_shortexp_dict_zip,timeobs_shortexp_dict_zip,FWHM_shortexp_dict_zip,mag_sex_matches_clean_shortexp_dict_zip,magerr_sex_matches_clean_shortexp_dict_zip,mag_atlas_matches_clean_shortexp_dict_zip,magerr_atlas_matches_clean_shortexp_dict_zip,mag_sex_matches_shortexp_dict_zip,magerr_sex_matches_shortexp_dict_zip,mag_atlas_matches_shortexp_dict_zip,magerr_atlas_matches_shortexp_dict_zip,zeropoints_all_shortexp_dict_zip,zeropoints_err_all_shortexp_dict_zip,airmasses_all_shortexp_dict_zip,timeobs_all_shortexp_dict_zip,FWHM_all_shortexp_dict_zip,donotstackscienceimages_shortexp_dict_zip=zip(*pool_mainzeropointcal_atlas.map(mainzeropointcal_mp,sciencefiles))

            logstring                                = collectziplists(logstring_zip)

            ######################## long exposure (120 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_longexp_dict          = collectzipdicts(stackscienceimages_longexp_dict_zip)
            zeropoints_longexp_dict                  = collectzipdicts(zeropoints_longexp_dict_zip)
            zeropoint_errs_longexp_dict              = collectzipdicts(zeropoint_errs_longexp_dict_zip)
            airmasses_longexp_dict                   = collectzipdicts(airmasses_longexp_dict_zip)
            timeobs_longexp_dict                     = collectzipdicts(timeobs_longexp_dict_zip)
            FWHM_longexp_dict                        = collectzipdicts(FWHM_longexp_dict_zip)
            # matched and cleaned
            mag_sex_matches_clean_longexp_dict       = collectzipdicts(mag_sex_matches_clean_longexp_dict_zip)
            magerr_sex_matches_clean_longexp_dict    = collectzipdicts(magerr_sex_matches_clean_longexp_dict_zip)
            mag_atlas_matches_clean_longexp_dict     = collectzipdicts(mag_atlas_matches_clean_longexp_dict_zip)
            magerr_atlas_matches_clean_longexp_dict  = collectzipdicts(magerr_atlas_matches_clean_longexp_dict_zip)
            # matched
            mag_sex_matches_longexp_dict             = collectzipdicts(mag_sex_matches_longexp_dict_zip)
            magerr_sex_matches_longexp_dict          = collectzipdicts(magerr_sex_matches_longexp_dict_zip)
            mag_atlas_matches_longexp_dict           = collectzipdicts(mag_atlas_matches_longexp_dict_zip)
            magerr_atlas_matches_longexp_dict        = collectzipdicts(magerr_atlas_matches_longexp_dict_zip)
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_longexp_dict              = collectzipdicts(zeropoints_all_longexp_dict_zip)
            zeropoints_err_all_longexp_dict          = collectzipdicts(zeropoints_err_all_longexp_dict_zip)
            airmasses_all_longexp_dict               = collectzipdicts(airmasses_all_longexp_dict_zip)
            timeobs_all_longexp_dict                 = collectzipdicts(timeobs_all_longexp_dict_zip)
            FWHM_all_longexp_dict                    = collectzipdicts(FWHM_all_longexp_dict_zip)
            donotstackscienceimages_longexp_dict     = collectzipdicts(donotstackscienceimages_longexp_dict_zip)

            ######################## short exposure (5 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_shortexp_dict         = collectzipdicts(stackscienceimages_shortexp_dict_zip)
            zeropoints_shortexp_dict                 = collectzipdicts(zeropoints_shortexp_dict_zip)
            zeropoint_errs_shortexp_dict             = collectzipdicts(zeropoint_errs_shortexp_dict_zip)
            airmasses_shortexp_dict                  = collectzipdicts(airmasses_shortexp_dict_zip)
            timeobs_shortexp_dict                    = collectzipdicts(timeobs_shortexp_dict_zip)
            FWHM_shortexp_dict                       = collectzipdicts(FWHM_shortexp_dict_zip)
            # matched and cleaned
            mag_sex_matches_clean_shortexp_dict      = collectzipdicts(mag_sex_matches_clean_shortexp_dict_zip)
            magerr_sex_matches_clean_shortexp_dict   = collectzipdicts(magerr_sex_matches_clean_shortexp_dict_zip)
            mag_atlas_matches_clean_shortexp_dict    = collectzipdicts(mag_atlas_matches_clean_shortexp_dict_zip)
            magerr_atlas_matches_clean_shortexp_dict = collectzipdicts(magerr_atlas_matches_clean_shortexp_dict_zip)
            # matched
            mag_sex_matches_shortexp_dict            = collectzipdicts(mag_sex_matches_shortexp_dict_zip)
            magerr_sex_matches_shortexp_dict         = collectzipdicts(magerr_sex_matches_shortexp_dict_zip)
            mag_atlas_matches_shortexp_dict          = collectzipdicts(mag_atlas_matches_shortexp_dict_zip)
            magerr_atlas_matches_shortexp_dict       = collectzipdicts(magerr_atlas_matches_shortexp_dict_zip)
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_shortexp_dict             = collectzipdicts(zeropoints_all_shortexp_dict_zip)
            zeropoints_err_all_shortexp_dict         = collectzipdicts(zeropoints_err_all_shortexp_dict_zip)
            airmasses_all_shortexp_dict              = collectzipdicts(airmasses_all_shortexp_dict_zip)
            timeobs_all_shortexp_dict                = collectzipdicts(timeobs_all_shortexp_dict_zip)
            FWHM_all_shortexp_dict                   = collectzipdicts(FWHM_all_shortexp_dict_zip)
            donotstackscienceimages_shortexp_dict    = collectzipdicts(donotstackscienceimages_shortexp_dict_zip)

        else:
            print "No multi-processing..."
            logstring                                = []

            ######################## long exposure (120 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_longexp_dict          = {}
            zeropoints_longexp_dict                  = {}
            zeropoint_errs_longexp_dict              = {}
            airmasses_longexp_dict                   = {}
            timeobs_longexp_dict                     = {}
            FWHM_longexp_dict                        = {}
            # matched and cleaned
            mag_sex_matches_clean_longexp_dict       = {}
            magerr_sex_matches_clean_longexp_dict    = {}
            mag_atlas_matches_clean_longexp_dict     = {}
            magerr_atlas_matches_clean_longexp_dict  = {}
            # matched
            mag_sex_matches_longexp_dict             = {}
            magerr_sex_matches_longexp_dict          = {}
            mag_atlas_matches_longexp_dict           = {}
            magerr_atlas_matches_longexp_dict        = {}
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_longexp_dict              = {}
            zeropoints_err_all_longexp_dict          = {}
            airmasses_all_longexp_dict               = {}
            timeobs_all_longexp_dict                 = {}
            FWHM_all_longexp_dict                    = {}
            donotstackscienceimages_longexp_dict     = {}

            ######################## short exposure (5 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_shortexp_dict         = {}
            zeropoints_shortexp_dict                 = {}
            zeropoint_errs_shortexp_dict             = {}
            airmasses_shortexp_dict                  = {}
            timeobs_shortexp_dict                    = {}
            FWHM_shortexp_dict                       = {}
            # matched and cleaned
            mag_sex_matches_clean_shortexp_dict      = {}
            magerr_sex_matches_clean_shortexp_dict   = {}
            mag_atlas_matches_clean_shortexp_dict    = {}
            magerr_atlas_matches_clean_shortexp_dict = {}
            # matched
            mag_sex_matches_shortexp_dict            = {}
            magerr_sex_matches_shortexp_dict         = {}
            mag_atlas_matches_shortexp_dict          = {}
            magerr_atlas_matches_shortexp_dict       = {}
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_shortexp_dict             = {}
            zeropoints_err_all_shortexp_dict         = {}
            airmasses_all_shortexp_dict              = {}
            timeobs_all_shortexp_dict                = {}
            FWHM_all_shortexp_dict                   = {}
            donotstackscienceimages_shortexp_dict    = {}

            i=0
            checkfiles = ["W3_1-S001-R001-C008-r_dupe-2-new-dsff.fts", "W3_1-S001-R001-C010-z_dupe-2-new-dsff.fts"]
            for file in sciencefiles:
                if file in checkfiles:
                    #if i<1:
                    logstring_temp,stackscienceimages_longexp_dict_temp,zeropoints_longexp_dict_temp,zeropoint_errs_longexp_dict_temp,airmasses_longexp_dict_temp,timeobs_longexp_dict_temp,FWHM_longexp_dict_temp,mag_sex_matches_clean_longexp_dict_temp,magerr_sex_matches_clean_longexp_dict_temp,mag_atlas_matches_clean_longexp_dict_temp,magerr_atlas_matches_clean_longexp_dict_temp,mag_sex_matches_longexp_dict_temp,magerr_sex_matches_longexp_dict_temp,mag_atlas_matches_longexp_dict_temp,magerr_atlas_matches_longexp_dict_temp,zeropoints_all_longexp_dict_temp,zeropoints_err_all_longexp_dict_temp,airmasses_all_longexp_dict_temp,timeobs_all_longexp_dict_temp,FWHM_all_longexp_dict_temp,donotstackscienceimages_longexp_dict_temp,stackscienceimages_shortexp_dict_temp,zeropoints_shortexp_dict_temp,zeropoint_errs_shortexp_dict_temp,airmasses_shortexp_dict_temp,timeobs_shortexp_dict_temp,FWHM_shortexp_dict_temp,mag_sex_matches_clean_shortexp_dict_temp,magerr_sex_matches_clean_shortexp_dict_temp,mag_atlas_matches_clean_shortexp_dict_temp,magerr_atlas_matches_clean_shortexp_dict_temp,mag_sex_matches_shortexp_dict_temp,magerr_sex_matches_shortexp_dict_temp,mag_atlas_matches_shortexp_dict_temp,magerr_atlas_matches_shortexp_dict_temp,zeropoints_all_shortexp_dict_temp,zeropoints_err_all_shortexp_dict_temp,airmasses_all_shortexp_dict_temp,timeobs_all_shortexp_dict_temp,FWHM_all_shortexp_dict_temp,donotstackscienceimages_shortexp_dict_temp=mainzeropointcal_atlas(atlasdata_trim,folder,folderdir,file,logbool=logbool,CHECKZP=CHECKZP)

                    logstring.append(logstring_temp)

                    ######################## long exposure (120 s) ########################
                    # zero slopes and 'good' zero point values
                    stackscienceimages_longexp_dict          = append_dict2dict(stackscienceimages_longexp_dict,stackscienceimages_longexp_dict_temp)
                    zeropoints_longexp_dict                  = append_dict2dict(zeropoints_longexp_dict,zeropoints_longexp_dict_temp)
                    zeropoint_errs_longexp_dict              = append_dict2dict(zeropoint_errs_longexp_dict,zeropoint_errs_longexp_dict_temp)
                    airmasses_longexp_dict                   = append_dict2dict(airmasses_longexp_dict,airmasses_longexp_dict_temp)
                    timeobs_longexp_dict                     = append_dict2dict(timeobs_longexp_dict,timeobs_longexp_dict_temp)
                    FWHM_longexp_dict                        = append_dict2dict(FWHM_longexp_dict,FWHM_longexp_dict_temp)
                    # matched and cleaned
                    mag_sex_matches_clean_longexp_dict       = append_dict2dict(mag_sex_matches_clean_longexp_dict,mag_sex_matches_clean_longexp_dict_temp)
                    magerr_sex_matches_clean_longexp_dict    = append_dict2dict(magerr_sex_matches_clean_longexp_dict,magerr_sex_matches_clean_longexp_dict_temp)
                    mag_atlas_matches_clean_longexp_dict     = append_dict2dict(mag_atlas_matches_clean_longexp_dict,mag_atlas_matches_clean_longexp_dict_temp)
                    magerr_atlas_matches_clean_longexp_dict  = append_dict2dict(magerr_atlas_matches_clean_longexp_dict,magerr_atlas_matches_clean_longexp_dict_temp)
                    # matched
                    mag_sex_matches_longexp_dict             = append_dict2dict(mag_sex_matches_longexp_dict,mag_sex_matches_longexp_dict_temp)
                    magerr_sex_matches_longexp_dict          = append_dict2dict(magerr_sex_matches_longexp_dict,magerr_sex_matches_longexp_dict_temp)
                    mag_atlas_matches_longexp_dict           = append_dict2dict(mag_atlas_matches_longexp_dict,mag_atlas_matches_longexp_dict_temp)
                    magerr_atlas_matches_longexp_dict        = append_dict2dict(magerr_atlas_matches_longexp_dict,magerr_atlas_matches_longexp_dict_temp)
                    # non-zero slopes and 'bad' zero-point values
                    zeropoints_all_longexp_dict              = append_dict2dict(zeropoints_all_longexp_dict,zeropoints_all_longexp_dict_temp)
                    zeropoints_err_all_longexp_dict          = append_dict2dict(zeropoints_err_all_longexp_dict,zeropoints_err_all_longexp_dict_temp)
                    airmasses_all_longexp_dict               = append_dict2dict(airmasses_all_longexp_dict,airmasses_all_longexp_dict_temp)
                    timeobs_all_longexp_dict                 = append_dict2dict(timeobs_all_longexp_dict,timeobs_all_longexp_dict_temp)
                    FWHM_all_longexp_dict                    = append_dict2dict(FWHM_all_longexp_dict,FWHM_all_longexp_dict_temp)
                    donotstackscienceimages_longexp_dict     = append_dict2dict(donotstackscienceimages_longexp_dict,donotstackscienceimages_longexp_dict_temp)

                    ######################## short exposure (5 s) ########################
                    # zero slopes and 'good' zero point values
                    stackscienceimages_shortexp_dict         = append_dict2dict(stackscienceimages_shortexp_dict,stackscienceimages_shortexp_dict_temp)
                    zeropoints_shortexp_dict                 = append_dict2dict(zeropoints_shortexp_dict,zeropoints_shortexp_dict_temp)
                    zeropoint_errs_shortexp_dict             = append_dict2dict(zeropoint_errs_shortexp_dict,zeropoint_errs_shortexp_dict_temp)
                    airmasses_shortexp_dict                  = append_dict2dict(airmasses_shortexp_dict,airmasses_shortexp_dict_temp)
                    timeobs_shortexp_dict                    = append_dict2dict(timeobs_shortexp_dict,timeobs_shortexp_dict_temp)
                    FWHM_shortexp_dict                       = append_dict2dict(FWHM_shortexp_dict,FWHM_shortexp_dict_temp)
                    # matched and cleaned
                    mag_sex_matches_clean_shortexp_dict      = append_dict2dict(mag_sex_matches_clean_shortexp_dict,mag_sex_matches_clean_shortexp_dict_temp)
                    magerr_sex_matches_clean_shortexp_dict   = append_dict2dict(magerr_sex_matches_clean_shortexp_dict,magerr_sex_matches_clean_shortexp_dict_temp)
                    mag_atlas_matches_clean_shortexp_dict    = append_dict2dict(mag_atlas_matches_clean_shortexp_dict,mag_atlas_matches_clean_shortexp_dict_temp)
                    magerr_atlas_matches_clean_shortexp_dict = append_dict2dict(magerr_atlas_matches_clean_shortexp_dict,magerr_atlas_matches_clean_shortexp_dict_temp)
                    # matched
                    mag_sex_matches_shortexp_dict            = append_dict2dict(mag_sex_matches_shortexp_dict,mag_sex_matches_shortexp_dict_temp)
                    magerr_sex_matches_shortexp_dict         = append_dict2dict(magerr_sex_matches_shortexp_dict,magerr_sex_matches_shortexp_dict_temp)
                    mag_atlas_matches_shortexp_dict          = append_dict2dict(mag_atlas_matches_shortexp_dict,mag_atlas_matches_shortexp_dict_temp)
                    magerr_atlas_matches_shortexp_dict       = append_dict2dict(magerr_atlas_matches_shortexp_dict,magerr_atlas_matches_shortexp_dict_temp)
                    # non-zero slopes and 'bad' zero-point values
                    zeropoints_all_shortexp_dict             = append_dict2dict(zeropoints_all_shortexp_dict,zeropoints_all_shortexp_dict_temp)
                    zeropoints_err_all_shortexp_dict         = append_dict2dict(zeropoints_err_all_shortexp_dict,zeropoints_err_all_shortexp_dict_temp)
                    airmasses_all_shortexp_dict              = append_dict2dict(airmasses_all_shortexp_dict,airmasses_all_shortexp_dict_temp)
                    timeobs_all_shortexp_dict                = append_dict2dict(timeobs_all_shortexp_dict,timeobs_all_shortexp_dict_temp)
                    FWHM_all_shortexp_dict                   = append_dict2dict(FWHM_all_shortexp_dict,FWHM_all_shortexp_dict_temp)
                    donotstackscienceimages_shortexp_dict    = append_dict2dict(donotstackscienceimages_shortexp_dict,donotstackscienceimages_shortexp_dict_temp)

                    i+=1

        ######################## long exposure (120 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_longexp_dict_flat         = flattendict(stackscienceimages_longexp_dict)
        print ">>>>>>>>>>>>>>> stackscienceimages_all_longexp_dict_flat (temp): <<<<<<<<<<<<<<"
        print stackscienceimages_all_longexp_dict_temp
        stackscienceimages_all_longexp_dict_temp     = append_dict2dict_temp(stackscienceimages_all_longexp_dict_temp,stackscienceimages_longexp_dict_flat)
        stackscienceimages_all_longexp_dict_flat     = flattendict(stackscienceimages_all_longexp_dict_temp,ARRAYS=True)
        print ">>>>>>>>>>>>>>> stackscienceimages_all_longexp_dict_flat (flat): <<<<<<<<<<<<<<"
        print stackscienceimages_all_longexp_dict_flat
        zeropoints_longexp_dict_flat                 = flattendict(zeropoints_longexp_dict)
        zeropoint_errs_longexp_dict_flat             = flattendict(zeropoint_errs_longexp_dict)
        airmasses_longexp_dict_flat                  = flattendict(airmasses_longexp_dict)
        timeobs_longexp_dict_flat                    = flattendict(timeobs_longexp_dict)
        FWHM_longexp_dict_flat                       = flattendict(FWHM_longexp_dict)
        # matched and cleaned
        mag_sex_matches_clean_longexp_dict_flat      = flattendict(mag_sex_matches_clean_longexp_dict,ARRAYS=True)
        magerr_sex_matches_clean_longexp_dict_flat   = flattendict(magerr_sex_matches_clean_longexp_dict,ARRAYS=True)
        mag_atlas_matches_clean_longexp_dict_flat    = flattendict(mag_atlas_matches_clean_longexp_dict,ARRAYS=True)
        magerr_atlas_matches_clean_longexp_dict_flat = flattendict(magerr_atlas_matches_clean_longexp_dict,ARRAYS=True)
        # matched
        mag_sex_matches_longexp_dict_flat            = flattendict(mag_sex_matches_longexp_dict,ARRAYS=True)
        magerr_sex_matches_longexp_dict_flat         = flattendict(magerr_sex_matches_longexp_dict,ARRAYS=True)
        mag_atlas_matches_longexp_dict_flat          = flattendict(mag_atlas_matches_longexp_dict,ARRAYS=True)
        magerr_atlas_matches_longexp_dict_flat       = flattendict(magerr_atlas_matches_longexp_dict,ARRAYS=True)
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_longexp_dict_flat             = flattendict(zeropoints_all_longexp_dict)
        zeropoints_err_all_longexp_dict_flat         = flattendict(zeropoints_err_all_longexp_dict)
        airmasses_all_longexp_dict_flat              = flattendict(airmasses_all_longexp_dict)
        timeobs_all_longexp_dict_flat                = flattendict(timeobs_all_longexp_dict)
        FWHM_all_longexp_dict_flat                   = flattendict(FWHM_all_longexp_dict)
        donotstackscienceimages_longexp_dict_flat    = flattendict(donotstackscienceimages_longexp_dict)

        ######################## short exposure (120 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_shortexp_dict_flat         = flattendict(stackscienceimages_shortexp_dict)
        print ">>>>>>>>>>>>>> stackscienceimages_all_shortexp_dict_flat (temp): <<<<<<<<<<<<<<"
        print stackscienceimages_all_shortexp_dict_temp
        stackscienceimages_all_shortexp_dict_temp     = append_dict2dict_temp(stackscienceimages_all_shortexp_dict_temp,stackscienceimages_shortexp_dict_flat)
        stackscienceimages_all_shortexp_dict_flat     = flattendict(stackscienceimages_all_shortexp_dict_temp,ARRAYS=True)
        print ">>>>>>>>>>>>>> stackscienceimages_all_shortexp_dict_flat (flat): <<<<<<<<<<<<<<"
        print stackscienceimages_all_shortexp_dict_flat
        zeropoints_shortexp_dict_flat                 = flattendict(zeropoints_shortexp_dict)
        zeropoint_errs_shortexp_dict_flat             = flattendict(zeropoint_errs_shortexp_dict)
        airmasses_shortexp_dict_flat                  = flattendict(airmasses_shortexp_dict)
        timeobs_shortexp_dict_flat                    = flattendict(timeobs_shortexp_dict)
        FWHM_shortexp_dict_flat                       = flattendict(FWHM_shortexp_dict)
        # matched and cleaned
        mag_sex_matches_clean_shortexp_dict_flat      = flattendict(mag_sex_matches_clean_shortexp_dict,ARRAYS=True)
        magerr_sex_matches_clean_shortexp_dict_flat   = flattendict(magerr_sex_matches_clean_shortexp_dict,ARRAYS=True)
        mag_atlas_matches_clean_shortexp_dict_flat    = flattendict(mag_atlas_matches_clean_shortexp_dict,ARRAYS=True)
        magerr_atlas_matches_clean_shortexp_dict_flat = flattendict(magerr_atlas_matches_clean_shortexp_dict,ARRAYS=True)
        # matched
        mag_sex_matches_shortexp_dict_flat            = flattendict(mag_sex_matches_shortexp_dict,ARRAYS=True)
        magerr_sex_matches_shortexp_dict_flat         = flattendict(magerr_sex_matches_shortexp_dict,ARRAYS=True)
        mag_atlas_matches_shortexp_dict_flat          = flattendict(mag_atlas_matches_shortexp_dict,ARRAYS=True)
        magerr_atlas_matches_shortexp_dict_flat       = flattendict(magerr_atlas_matches_shortexp_dict,ARRAYS=True)
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_shortexp_dict_flat             = flattendict(zeropoints_all_shortexp_dict)
        zeropoints_err_all_shortexp_dict_flat         = flattendict(zeropoints_err_all_shortexp_dict)
        airmasses_all_shortexp_dict_flat              = flattendict(airmasses_all_shortexp_dict)
        timeobs_all_shortexp_dict_flat                = flattendict(timeobs_all_shortexp_dict)
        FWHM_all_shortexp_dict_flat                   = flattendict(FWHM_all_shortexp_dict)
        donotstackscienceimages_shortexp_dict_flat    = flattendict(donotstackscienceimages_shortexp_dict)

        if VERBOSE:
            print "Plotting all zero-point data for each passband..."

        # plot zeropoint data for each passband
        plotzp_all(mag_sex_matches_clean_longexp_dict_flat,magerr_sex_matches_clean_longexp_dict_flat,mag_atlas_matches_clean_longexp_dict_flat,magerr_atlas_matches_clean_longexp_dict_flat,mag_sex_matches_longexp_dict_flat,magerr_sex_matches_longexp_dict_flat,mag_atlas_matches_longexp_dict_flat,magerr_atlas_matches_longexp_dict_flat,atlasplotfolderdir,"LONG",MULTIPROCESSING_CAL,CHECKZP)
        plotzp_all(mag_sex_matches_clean_shortexp_dict_flat,magerr_sex_matches_clean_shortexp_dict_flat,mag_atlas_matches_clean_shortexp_dict_flat,magerr_atlas_matches_clean_shortexp_dict_flat,mag_sex_matches_shortexp_dict_flat,magerr_sex_matches_shortexp_dict_flat,mag_atlas_matches_shortexp_dict_flat,magerr_atlas_matches_shortexp_dict_flat,atlasplotfolderdir,"SHORT",MULTIPROCESSING_CAL,CHECKZP)        

        if VERBOSE:
            print "Plotting zero-point correlations for each passband..."

        # plot zeropoint correlations
        plotzpcorrs(zeropoints_longexp_dict_flat,zeropoint_errs_longexp_dict_flat,zeropoints_all_longexp_dict_flat,zeropoints_err_all_longexp_dict_flat,airmasses_longexp_dict_flat,airmasses_all_longexp_dict_flat,timeobs_longexp_dict_flat,timeobs_all_longexp_dict_flat,FWHM_longexp_dict_flat,FWHM_all_longexp_dict_flat,"LONG",atlasplotfolderdir,CHECKZP)
        plotzpcorrs(zeropoints_shortexp_dict_flat,zeropoint_errs_shortexp_dict_flat,zeropoints_all_shortexp_dict_flat,zeropoints_err_all_shortexp_dict_flat,airmasses_shortexp_dict_flat,airmasses_all_shortexp_dict_flat,timeobs_shortexp_dict_flat,timeobs_all_shortexp_dict_flat,FWHM_shortexp_dict_flat,FWHM_all_shortexp_dict_flat,"SHORT",atlasplotfolderdir,CHECKZP)

        # clear variables from memory (running the initial calibration and the second calibration step back-to-back will double their memory and raise a memory error)
        if VERBOSE:
            print "clearing calibration vairables to minimize memory..."

        ######################## long exposure (120 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_longexp_dict = {}             # non-flat dictionaries
        zeropoints_longexp_dict = {}                     # non-flat dictionaries
        zeropoint_errs_longexp_dict = {}                 # non-flat dictionaries
        airmasses_longexp_dict_flat = {}                 # for ZP correlation plots
        airmasses_longexp_dict = {}                      # for ZP correlation plots
        timeobs_longexp_dict_flat = {}                   # for ZP correlation plots
        timeobs_longexp_dict = {}                        # for ZP correlation plots
        FWHM_longexp_dict_flat = {}                      # for ZP correlation plots
        FWHM_longexp_dict = {}                           # for ZP correlation plots
        # matched and cleaned
        mag_sex_matches_clean_longexp_dict_flat = {}     # for ZP plots
        mag_sex_matches_clean_longexp_dict = {}          # for ZP plots
        magerr_sex_matches_clean_longexp_dict_flat = {}  # for ZP plots
        magerr_sex_matches_clean_longexp_dict = {}       # for ZP plots
        mag_ps_matches_clean_longexp_dict_flat = {}      # for ZP plots
        mag_ps_matches_clean_longexp_dict = {}           # for ZP plots
        magerr_ps_matches_clean_longexp_dict_flat = {}   # for ZP plots
        magerr_ps_matches_clean_longexp_dict = {}        # for ZP plots
        # matched
        mag_sex_matches_longexp_dict_flat = {}           # for ZP plots
        mag_sex_matches_longexp_dict = {}                # for ZP plots
        magerr_sex_matches_longexp_dict_flat = {}        # for ZP plots
        magerr_sex_matches_longexp_dict = {}             # for ZP plots
        mag_ps_matches_longexp_dict_flat = {}            # for ZP plots
        mag_ps_matches_longexp_dict = {}                 # for ZP plots
        magerr_ps_matches_longexp_dict_flat = {}         # for ZP plots
        magerr_ps_matches_longexp_dict = {}              # for ZP plots
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_longexp_dict_flat = {}            # for ZP plots
        zeropoints_all_longexp_dict = {}                 # for ZP plots
        zeropoints_err_all_longexp_dict_flat = {}        # for ZP plots
        zeropoints_err_all_longexp_dict = {}             # for ZP plots
        airmasses_all_longexp_dict_flat = {}             # for ZP correlation plots
        airmasses_all_longexp_dict = {}                  # for ZP correlation plots
        timeobs_all_longexp_dict_flat = {}               # for ZP correlation plots
        timeobs_all_longexp_dict = {}                    # for ZP correlation plots
        FWHM_all_longexp_dict_flat = {}                  # for ZP correlation plots
        FWHM_all_longexp_dict = {}                       # for ZP correlation plots
        donotstackscienceimages_longexp_dict = {}        # non-flat dictionary
        
        ######################## short exposure (5 s) ########################
        # zero slopes and 'good' zero point values
        stackscienceimages_shortexp_dict = {}            # non-flat dictionaries
        zeropoints_shortexp_dict = {}                    # non-flat dictionaries
        zeropoint_errs_shortexp_dict = {}                # non-flat dictionaries
        airmasses_shortexp_dict_flat = {}                # for ZP correlation plots
        airmasses_shortexp_dict = {}                     # for ZP correlation plots
        timeobs_shortexp_dict_flat = {}                  # for ZP correlation plots
        timeobs_shortexp_dict = {}                       # for ZP correlation plots
        FWHM_shortexp_dict_flat = {}                     # for ZP correlation plots
        FWHM_shortexp_dict = {}                          # for ZP correlation plots
        # matched and cleaned
        mag_sex_matches_clean_shortexp_dict_flat = {}    # for ZP plots
        mag_sex_matches_clean_shortexp_dict = {}         # for ZP plots
        magerr_sex_matches_clean_shortexp_dict_flat = {} # for ZP plots
        magerr_sex_matches_clean_shortexp_dict = {}      # for ZP plots
        mag_ps_matches_clean_shortexp_dict_flat = {}     # for ZP plots
        mag_ps_matches_clean_shortexp_dict = {}          # for ZP plots
        magerr_ps_matches_clean_shortexp_dict_flat = {}  # for ZP plots
        magerr_ps_matches_clean_shortexp_dict = {}       # for ZP plots
        # matched
        mag_sex_matches_shortexp_dict_flat = {}          # for ZP plots
        mag_sex_matches_shortexp_dict = {}               # for ZP plots
        magerr_sex_matches_shortexp_dict_flat = {}       # for ZP plots
        magerr_sex_matches_shortexp_dict = {}            # for ZP plots
        mag_ps_matches_shortexp_dict_flat = {}           # for ZP plots
        mag_ps_matches_shortexp_dict = {}                # for ZP plots
        magerr_ps_matches_shortexp_dict_flat = {}        # for ZP plots
        magerr_ps_matches_shortexp_dict = {}             # for ZP plots
        # non-zero slopes and 'bad' zero-point values
        zeropoints_all_shortexp_dict_flat = {}           # for ZP plots
        zeropoints_all_shortexp_dict = {}                # for ZP plots
        zeropoints_err_all_shortexp_dict_flat = {}       # for ZP plots
        zeropoints_err_all_shortexp_dict = {}            # for ZP plots
        airmasses_all_shortexp_dict_flat = {}            # for ZP correlation plots
        airmasses_all_shortexp_dict = {}                 # for ZP correlation plots
        timeobs_all_shortexp_dict_flat = {}              # for ZP correlation plots
        timeobs_all_shortexp_dict = {}                   # for ZP correlation plots
        FWHM_all_shortexp_dict_flat = {}                 # for ZP correlation plots
        FWHM_all_shortexp_dict = {}                      # for ZP correlation plots
        donotstackscienceimages_shortexp_dict = {}       # non-flat dictionary

        print "\nscience images to be stacked (long exp):"
        print stackscienceimages_longexp_dict_flat
        print "\nzero points (long exp):"
        print zeropoints_longexp_dict
        print "\nzero point errors (long exp):"
        print zeropoint_errs_longexp_dict
        print "\nbad science images (long exposure):"
        print donotstackscienceimages_longexp_dict_flat
        print "\nscience images to be stacked (short exp):"
        print stackscienceimages_shortexp_dict_flat
        print "\nbad science images (short exposure):"
        print donotstackscienceimages_shortexp_dict_flat

        if LOGOUTPUT:
            # write calibration log
            if VERBOSE:
                print "writing zeropoint log"
            f=open(logdir+"zeropoint.log","w+")
            for item in logstring:
                f.write("%s\n" % item)
            f.close()

    print ">>>>>>>>>>>>>> stackscienceimages_all_longexp_dict_flat (final) <<<<<<<<<<<<<<<<"
    print stackscienceimages_all_longexp_dict_flat
    print ">>>>>>>>>>>>>> stackscienceimages_all_shortexp_dict_flat (final) <<<<<<<<<<<<<<<<"
    print stackscienceimages_all_shortexp_dict_flat
    # store dictionary of files to stack
    fwritedict_pickle(stackscienceimages_all_longexp_dict_flat,"stackscienceimages_all_longexp_dict_pickle.txt")
    fwritedict_pickle(stackscienceimages_all_shortexp_dict_flat,"stackscienceimages_all_shortexp_dict_pickle.txt")
    fwritedict_text(stackscienceimages_all_longexp_dict_flat,"stackscienceimages_all_longexp_dict.txt")
    fwritedict_text(stackscienceimages_all_shortexp_dict_flat,"stackscienceimages_all_shortexp_dict.txt")

else:
    if VERBOSE:
        print "\n---------------------------- skipping ATLAS calilbration -------------------------"


#########################################################################################################
#                                                                                                       #
#                                  CHECK ZEROPOINT CALIBRATION                                          #
#                                                                                                       #
#########################################################################################################

#########################################################################################################
#                                                                                                       #
#                               CHECK ATLAS ZEROPOINT CALIBRATION                                       #
#                                                                                                       #
#########################################################################################################

if ATLAS and CHECKCAL:
    CHECKZP=True

    # read in catalogue data
    atlascatf = atlasdir+atlascat

    if VERBOSE:
        print "\n--------------------------- checking ATLAS calilbration --------------------------"
    print "Reading in ATLAS catalogue..."
    # trimmed version of ATLAS data
    atlastrimcatf   = atlasdir+atlastrimcat
    hdul_atlastrim  = fits.open(atlastrimcatf)
    atlasdata_trim  = hdul_atlastrim[1].data

    #################################### FIND MAGNITUDE ZERO POINT ###########################################

    atlascalibrateddir       = calibrateddir+"atlas/"
    atlasplotdir_calibrated  = calibratedplotdir+"atlas/"

    if pathexists(atlasplotdir_calibrated):
        os.system("rm -rf "+atlasplotdir_calibrated)
    if VERBOSE:
        print "cleared directory "+atlasplotdir_calibrated

    os.system("mkdir "+atlasplotdir_calibrated)
    if VERBOSE:
        print "created directory "+atlasplotdir_calibrated

    for folder in folders:
        atlascalibratedplotfolderdir = atlasplotdir_calibrated+folder+"/"
        if pathexists(atlascalibratedplotfolderdir):
            os.system("rm -rf "+atlascalibratedplotfolderdir)
            if VERBOSE:
                print "cleared directory "+atlascalibratedplotfolderdir
        os.system("mkdir "+atlascalibratedplotfolderdir)
        if VERBOSE:
            print "created directory "+atlascalibratedplotfolderdir

    if SEXTRACTOR:
        checkdirs = [sexobjcheckzpdir,sexcatcheckzpdir,sexbgcheckzpdir]

        # clear catalog directories
        if pathexists(sexcheckzpdir):
            os.system("rm -rf "+sexcheckzpdir)
            if VERBOSE:
                print "cleared directory "+sexcheckzpdir
        os.system("mkdir "+sexcheckzpdir)
        for dir in checkdirs:
            os.system("mkdir "+dir)
            for folder in folders:
                os.system("mkdir "+dir+folder)

    if LOGOUTPUT:
        logbool=True
        logstring=[]
    else:
        logbool=False

    pool_mainzeropointcal_atlas = multiprocessing.Pool(ncores_atlas)

    if VERBOSE:
        print "\n----------------- obtaining overlap between ATLAS and SExtractor -----------------"
    if LOGOUTPUT:
        logstring.append("----------------- obtaining overlap between ATLAS and SExtractor -----------------")

    for folder in folders:
        calibratedfolderdir          = atlascalibrateddir+folder+"/"
        atlascalibratedplotfolderdir = atlasplotdir_calibrated+folder+"/"
        if VERBOSE:
            print "processing folder "+folder+"..."
        if LOGOUTPUT:
            logstring.append(" > processing folder "+str(folder)+"...")
        #folderdir = findfolderdir(folder)
        if VERBOSE:
            print "using data directory "+calibratedfolderdir+" for folder "+folder
        if LOGOUTPUT:
            logstring.append("using data directory "+str(calibratedfolderdir)+" for folder "+str(folder))
        sciencefiles = os.listdir(calibratedfolderdir)

        # only proceed if there if usable data
        if len(sciencefiles)>0:
            if VERBOSE:
                print "Using science files:", sciencefiles
            if MULTIPROCESSING_CAL:
                print "Multiprocessing..."
                mainzeropointcal_mp=partial(mainzeropointcal_atlas,atlasdata_trim,folder,calibratedfolderdir,logbool=logbool,CHECKZP=CHECKZP)
                logstring_zip,stackscienceimages_longexp_dict_zip,zeropoints_longexp_dict_zip,zeropoint_errs_longexp_dict_zip,airmasses_longexp_dict_zip,timeobs_longexp_dict_zip,FWHM_longexp_dict_zip,mag_sex_matches_clean_longexp_dict_zip,magerr_sex_matches_clean_longexp_dict_zip,mag_atlas_matches_clean_longexp_dict_zip,magerr_atlas_matches_clean_longexp_dict_zip,mag_sex_matches_longexp_dict_zip,magerr_sex_matches_longexp_dict_zip,mag_atlas_matches_longexp_dict_zip,magerr_atlas_matches_longexp_dict_zip,zeropoints_all_longexp_dict_zip,zeropoints_err_all_longexp_dict_zip,airmasses_all_longexp_dict_zip,timeobs_all_longexp_dict_zip,FWHM_all_longexp_dict_zip,donotstackscienceimages_longexp_dict_zip,stackscienceimages_shortexp_dict_zip,zeropoints_shortexp_dict_zip,zeropoint_errs_shortexp_dict_zip,airmasses_shortexp_dict_zip,timeobs_shortexp_dict_zip,FWHM_shortexp_dict_zip,mag_sex_matches_clean_shortexp_dict_zip,magerr_sex_matches_clean_shortexp_dict_zip,mag_atlas_matches_clean_shortexp_dict_zip,magerr_atlas_matches_clean_shortexp_dict_zip,mag_sex_matches_shortexp_dict_zip,magerr_sex_matches_shortexp_dict_zip,mag_atlas_matches_shortexp_dict_zip,magerr_atlas_matches_shortexp_dict_zip,zeropoints_all_shortexp_dict_zip,zeropoints_err_all_shortexp_dict_zip,airmasses_all_shortexp_dict_zip,timeobs_all_shortexp_dict_zip,FWHM_all_shortexp_dict_zip,donotstackscienceimages_shortexp_dict_zip=zip(*pool_mainzeropointcal_atlas.map(mainzeropointcal_mp,sciencefiles))

                logstring                                = collectziplists(logstring_zip)

                ######################## long exposure (120 s) ########################
                # zero slopes and 'good' zero point values
                stackscienceimages_longexp_dict          = collectzipdicts(stackscienceimages_longexp_dict_zip)
                zeropoints_longexp_dict                  = collectzipdicts(zeropoints_longexp_dict_zip)
                zeropoint_errs_longexp_dict              = collectzipdicts(zeropoint_errs_longexp_dict_zip)
                airmasses_longexp_dict                   = collectzipdicts(airmasses_longexp_dict_zip)
                timeobs_longexp_dict                     = collectzipdicts(timeobs_longexp_dict_zip)
                FWHM_longexp_dict                        = collectzipdicts(FWHM_longexp_dict_zip)
                # matched and cleaned
                mag_sex_matches_clean_longexp_dict       = collectzipdicts(mag_sex_matches_clean_longexp_dict_zip)
                magerr_sex_matches_clean_longexp_dict    = collectzipdicts(magerr_sex_matches_clean_longexp_dict_zip)
                mag_atlas_matches_clean_longexp_dict     = collectzipdicts(mag_atlas_matches_clean_longexp_dict_zip)
                magerr_atlas_matches_clean_longexp_dict  = collectzipdicts(magerr_atlas_matches_clean_longexp_dict_zip)
                # matched
                mag_sex_matches_longexp_dict             = collectzipdicts(mag_sex_matches_longexp_dict_zip)
                magerr_sex_matches_longexp_dict          = collectzipdicts(magerr_sex_matches_longexp_dict_zip)
                mag_atlas_matches_longexp_dict           = collectzipdicts(mag_atlas_matches_longexp_dict_zip)
                magerr_atlas_matches_longexp_dict        = collectzipdicts(magerr_atlas_matches_longexp_dict_zip)
                # non-zero slopes and 'bad' zero-point values
                zeropoints_all_longexp_dict              = collectzipdicts(zeropoints_all_longexp_dict_zip)
                zeropoints_err_all_longexp_dict          = collectzipdicts(zeropoints_err_all_longexp_dict_zip)
                airmasses_all_longexp_dict               = collectzipdicts(airmasses_all_longexp_dict_zip)
                timeobs_all_longexp_dict                 = collectzipdicts(timeobs_all_longexp_dict_zip)
                FWHM_all_longexp_dict                    = collectzipdicts(FWHM_all_longexp_dict_zip)
                donotstackscienceimages_longexp_dict     = collectzipdicts(donotstackscienceimages_longexp_dict_zip)

                ######################## short exposure (5 s) ########################
                # zero slopes and 'good' zero point values
                stackscienceimages_shortexp_dict         = collectzipdicts(stackscienceimages_shortexp_dict_zip)
                zeropoints_shortexp_dict                 = collectzipdicts(zeropoints_shortexp_dict_zip)
                zeropoint_errs_shortexp_dict             = collectzipdicts(zeropoint_errs_shortexp_dict_zip)
                airmasses_shortexp_dict                  = collectzipdicts(airmasses_shortexp_dict_zip)
                timeobs_shortexp_dict                    = collectzipdicts(timeobs_shortexp_dict_zip)
                FWHM_shortexp_dict                       = collectzipdicts(FWHM_shortexp_dict_zip)
                # matched and cleaned
                mag_sex_matches_clean_shortexp_dict      = collectzipdicts(mag_sex_matches_clean_shortexp_dict_zip)
                magerr_sex_matches_clean_shortexp_dict   = collectzipdicts(magerr_sex_matches_clean_shortexp_dict_zip)
                mag_atlas_matches_clean_shortexp_dict    = collectzipdicts(mag_atlas_matches_clean_shortexp_dict_zip)
                magerr_atlas_matches_clean_shortexp_dict = collectzipdicts(magerr_atlas_matches_clean_shortexp_dict_zip)
                # matched
                mag_sex_matches_shortexp_dict            = collectzipdicts(mag_sex_matches_shortexp_dict_zip)
                magerr_sex_matches_shortexp_dict         = collectzipdicts(magerr_sex_matches_shortexp_dict_zip)
                mag_atlas_matches_shortexp_dict          = collectzipdicts(mag_atlas_matches_shortexp_dict_zip)
                magerr_atlas_matches_shortexp_dict       = collectzipdicts(magerr_atlas_matches_shortexp_dict_zip)
                # non-zero slopes and 'bad' zero-point values
                zeropoints_all_shortexp_dict             = collectzipdicts(zeropoints_all_shortexp_dict_zip)
                zeropoints_err_all_shortexp_dict         = collectzipdicts(zeropoints_err_all_shortexp_dict_zip)
                airmasses_all_shortexp_dict              = collectzipdicts(airmasses_all_shortexp_dict_zip)
                timeobs_all_shortexp_dict                = collectzipdicts(timeobs_all_shortexp_dict_zip)
                FWHM_all_shortexp_dict                   = collectzipdicts(FWHM_all_shortexp_dict_zip)
                donotstackscienceimages_shortexp_dict    = collectzipdicts(donotstackscienceimages_shortexp_dict_zip)

            else:
                print "No multi-processing..."
                logstring                                = []

                ######################## long exposure (120 s) ########################
                # zero slopes and 'good' zero point values
                stackscienceimages_longexp_dict          = {}
                zeropoints_longexp_dict                  = {}
                zeropoint_errs_longexp_dict              = {}
                airmasses_longexp_dict                   = {}
                timeobs_longexp_dict                     = {}
                FWHM_longexp_dict                        = {}
                # matched and cleaned
                mag_sex_matches_clean_longexp_dict       = {}
                magerr_sex_matches_clean_longexp_dict    = {}
                mag_atlas_matches_clean_longexp_dict     = {}
                magerr_atlas_matches_clean_longexp_dict  = {}
                # matched
                mag_sex_matches_longexp_dict             = {}
                magerr_sex_matches_longexp_dict          = {}
                mag_atlas_matches_longexp_dict           = {}
                magerr_atlas_matches_longexp_dict        = {}
                # non-zero slopes and 'bad' zero-point values
                zeropoints_all_longexp_dict              = {}
                zeropoints_err_all_longexp_dict          = {}
                airmasses_all_longexp_dict               = {}
                timeobs_all_longexp_dict                 = {}
                FWHM_all_longexp_dict                    = {}
                donotstackscienceimages_longexp_dict     = {}

                ######################## short exposure (5 s) ########################
                # zero slopes and 'good' zero point values
                stackscienceimages_shortexp_dict         = {}
                zeropoints_shortexp_dict                 = {}
                zeropoint_errs_shortexp_dict             = {}
                airmasses_shortexp_dict                  = {}
                timeobs_shortexp_dict                    = {}
                FWHM_shortexp_dict                       = {}
                # matched and cleaned
                mag_sex_matches_clean_shortexp_dict      = {}
                magerr_sex_matches_clean_shortexp_dict   = {}
                mag_atlas_matches_clean_shortexp_dict    = {}
                magerr_atlas_matches_clean_shortexp_dict = {}
                # matched
                mag_sex_matches_shortexp_dict            = {}
                magerr_sex_matches_shortexp_dict         = {}
                mag_atlas_matches_shortexp_dict          = {}
                magerr_atlas_matches_shortexp_dict       = {}
                # non-zero slopes and 'bad' zero-point values
                zeropoints_all_shortexp_dict             = {}
                zeropoints_err_all_shortexp_dict         = {}
                airmasses_all_shortexp_dict              = {}
                timeobs_all_shortexp_dict                = {}
                FWHM_all_shortexp_dict                   = {}
                donotstackscienceimages_shortexp_dict    = {}

                i=0
                checkfiles = ["W3_1-S001-R001-C008-r_dupe-2-new-dsff_calibrated.fts", "W3_1-S001-R001-C010-z_dupe-2-new-dsff_calibrated.fts"]
                for file in sciencefiles:
                    if file in checkfiles:
                        #if i<1:
                        logstring_temp,stackscienceimages_longexp_dict_temp,zeropoints_longexp_dict_temp,zeropoint_errs_longexp_dict_temp,airmasses_longexp_dict_temp,timeobs_longexp_dict_temp,FWHM_longexp_dict_temp,mag_sex_matches_clean_longexp_dict_temp,magerr_sex_matches_clean_longexp_dict_temp,mag_atlas_matches_clean_longexp_dict_temp,magerr_atlas_matches_clean_longexp_dict_temp,mag_sex_matches_longexp_dict_temp,magerr_sex_matches_longexp_dict_temp,mag_atlas_matches_longexp_dict_temp,magerr_atlas_matches_longexp_dict_temp,zeropoints_all_longexp_dict_temp,zeropoints_err_all_longexp_dict_temp,airmasses_all_longexp_dict_temp,timeobs_all_longexp_dict_temp,FWHM_all_longexp_dict_temp,donotstackscienceimages_longexp_dict_temp,stackscienceimages_shortexp_dict_temp,zeropoints_shortexp_dict_temp,zeropoint_errs_shortexp_dict_temp,airmasses_shortexp_dict_temp,timeobs_shortexp_dict_temp,FWHM_shortexp_dict_temp,mag_sex_matches_clean_shortexp_dict_temp,magerr_sex_matches_clean_shortexp_dict_temp,mag_atlas_matches_clean_shortexp_dict_temp,magerr_atlas_matches_clean_shortexp_dict_temp,mag_sex_matches_shortexp_dict_temp,magerr_sex_matches_shortexp_dict_temp,mag_atlas_matches_shortexp_dict_temp,magerr_atlas_matches_shortexp_dict_temp,zeropoints_all_shortexp_dict_temp,zeropoints_err_all_shortexp_dict_temp,airmasses_all_shortexp_dict_temp,timeobs_all_shortexp_dict_temp,FWHM_all_shortexp_dict_temp,donotstackscienceimages_shortexp_dict_temp=mainzeropointcal_atlas(atlasdata_trim,folder,calibratedfolderdir,file,logbool=logbool,CHECKZP=CHECKZP)

                        logstring.append(logstring_temp)

                        ######################## long exposure (120 s) ########################
                        # zero slopes and 'good' zero point values
                        stackscienceimages_longexp_dict          = append_dict2dict(stackscienceimages_longexp_dict,stackscienceimages_longexp_dict_temp)
                        zeropoints_longexp_dict                  = append_dict2dict(zeropoints_longexp_dict,zeropoints_longexp_dict_temp)
                        zeropoint_errs_longexp_dict              = append_dict2dict(zeropoint_errs_longexp_dict,zeropoint_errs_longexp_dict_temp)
                        airmasses_longexp_dict                   = append_dict2dict(airmasses_longexp_dict,airmasses_longexp_dict_temp)
                        timeobs_longexp_dict                     = append_dict2dict(timeobs_longexp_dict,timeobs_longexp_dict_temp)
                        FWHM_longexp_dict                        = append_dict2dict(FWHM_longexp_dict,FWHM_longexp_dict_temp)
                        # matched and cleaned
                        mag_sex_matches_clean_longexp_dict       = append_dict2dict(mag_sex_matches_clean_longexp_dict,mag_sex_matches_clean_longexp_dict_temp)
                        magerr_sex_matches_clean_longexp_dict    = append_dict2dict(magerr_sex_matches_clean_longexp_dict,magerr_sex_matches_clean_longexp_dict_temp)
                        mag_atlas_matches_clean_longexp_dict     = append_dict2dict(mag_atlas_matches_clean_longexp_dict,mag_atlas_matches_clean_longexp_dict_temp)
                        magerr_atlas_matches_clean_longexp_dict  = append_dict2dict(magerr_atlas_matches_clean_longexp_dict,magerr_atlas_matches_clean_longexp_dict_temp)
                        # matched
                        mag_sex_matches_longexp_dict             = append_dict2dict(mag_sex_matches_longexp_dict,mag_sex_matches_longexp_dict_temp)
                        magerr_sex_matches_longexp_dict          = append_dict2dict(magerr_sex_matches_longexp_dict,magerr_sex_matches_longexp_dict_temp)
                        mag_atlas_matches_longexp_dict           = append_dict2dict(mag_atlas_matches_longexp_dict,mag_atlas_matches_longexp_dict_temp)
                        magerr_atlas_matches_longexp_dict        = append_dict2dict(magerr_atlas_matches_longexp_dict,magerr_atlas_matches_longexp_dict_temp)
                        # non-zero slopes and 'bad' zero-point values
                        zeropoints_all_longexp_dict              = append_dict2dict(zeropoints_all_longexp_dict,zeropoints_all_longexp_dict_temp)
                        zeropoints_err_all_longexp_dict          = append_dict2dict(zeropoints_err_all_longexp_dict,zeropoints_err_all_longexp_dict_temp)
                        airmasses_all_longexp_dict               = append_dict2dict(airmasses_all_longexp_dict,airmasses_all_longexp_dict_temp)
                        timeobs_all_longexp_dict                 = append_dict2dict(timeobs_all_longexp_dict,timeobs_all_longexp_dict_temp)
                        FWHM_all_longexp_dict                    = append_dict2dict(FWHM_all_longexp_dict,FWHM_all_longexp_dict_temp)
                        donotstackscienceimages_longexp_dict     = append_dict2dict(donotstackscienceimages_longexp_dict,donotstackscienceimages_longexp_dict_temp)

                        ######################## short exposure (5 s) ########################
                        # zero slopes and 'good' zero point values
                        stackscienceimages_shortexp_dict         = append_dict2dict(stackscienceimages_shortexp_dict,stackscienceimages_shortexp_dict_temp)
                        zeropoints_shortexp_dict                 = append_dict2dict(zeropoints_shortexp_dict,zeropoints_shortexp_dict_temp)
                        zeropoint_errs_shortexp_dict             = append_dict2dict(zeropoint_errs_shortexp_dict,zeropoint_errs_shortexp_dict_temp)
                        airmasses_shortexp_dict                  = append_dict2dict(airmasses_shortexp_dict,airmasses_shortexp_dict_temp)
                        timeobs_shortexp_dict                    = append_dict2dict(timeobs_shortexp_dict,timeobs_shortexp_dict_temp)
                        FWHM_shortexp_dict                       = append_dict2dict(FWHM_shortexp_dict,FWHM_shortexp_dict_temp)
                        # matched and cleaned
                        mag_sex_matches_clean_shortexp_dict      = append_dict2dict(mag_sex_matches_clean_shortexp_dict,mag_sex_matches_clean_shortexp_dict_temp)
                        magerr_sex_matches_clean_shortexp_dict   = append_dict2dict(magerr_sex_matches_clean_shortexp_dict,magerr_sex_matches_clean_shortexp_dict_temp)
                        mag_atlas_matches_clean_shortexp_dict    = append_dict2dict(mag_atlas_matches_clean_shortexp_dict,mag_atlas_matches_clean_shortexp_dict_temp)
                        magerr_atlas_matches_clean_shortexp_dict = append_dict2dict(magerr_atlas_matches_clean_shortexp_dict,magerr_atlas_matches_clean_shortexp_dict_temp)
                        # matched
                        mag_sex_matches_shortexp_dict            = append_dict2dict(mag_sex_matches_shortexp_dict,mag_sex_matches_shortexp_dict_temp)
                        magerr_sex_matches_shortexp_dict         = append_dict2dict(magerr_sex_matches_shortexp_dict,magerr_sex_matches_shortexp_dict_temp)
                        mag_atlas_matches_shortexp_dict          = append_dict2dict(mag_atlas_matches_shortexp_dict,mag_atlas_matches_shortexp_dict_temp)
                        magerr_atlas_matches_shortexp_dict       = append_dict2dict(magerr_atlas_matches_shortexp_dict,magerr_atlas_matches_shortexp_dict_temp)
                        # non-zero slopes and 'bad' zero-point values
                        zeropoints_all_shortexp_dict             = append_dict2dict(zeropoints_all_shortexp_dict,zeropoints_all_shortexp_dict_temp)
                        zeropoints_err_all_shortexp_dict         = append_dict2dict(zeropoints_err_all_shortexp_dict,zeropoints_err_all_shortexp_dict_temp)
                        airmasses_all_shortexp_dict              = append_dict2dict(airmasses_all_shortexp_dict,airmasses_all_shortexp_dict_temp)
                        timeobs_all_shortexp_dict                = append_dict2dict(timeobs_all_shortexp_dict,timeobs_all_shortexp_dict_temp)
                        FWHM_all_shortexp_dict                   = append_dict2dict(FWHM_all_shortexp_dict,FWHM_all_shortexp_dict_temp)
                        donotstackscienceimages_shortexp_dict    = append_dict2dict(donotstackscienceimages_shortexp_dict,donotstackscienceimages_shortexp_dict_temp)

                        i+=1

            ######################## long exposure (120 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_longexp_dict_flat         = flattendict(stackscienceimages_longexp_dict)
            zeropoints_longexp_dict_flat                 = flattendict(zeropoints_longexp_dict)
            zeropoint_errs_longexp_dict_flat             = flattendict(zeropoint_errs_longexp_dict)
            airmasses_longexp_dict_flat                  = flattendict(airmasses_longexp_dict)
            timeobs_longexp_dict_flat                    = flattendict(timeobs_longexp_dict)
            FWHM_longexp_dict_flat                       = flattendict(FWHM_longexp_dict)
            # matched and cleaned
            mag_sex_matches_clean_longexp_dict_flat      = flattendict(mag_sex_matches_clean_longexp_dict,ARRAYS=True)
            magerr_sex_matches_clean_longexp_dict_flat   = flattendict(magerr_sex_matches_clean_longexp_dict,ARRAYS=True)
            mag_atlas_matches_clean_longexp_dict_flat    = flattendict(mag_atlas_matches_clean_longexp_dict,ARRAYS=True)
            magerr_atlas_matches_clean_longexp_dict_flat = flattendict(magerr_atlas_matches_clean_longexp_dict,ARRAYS=True)
            # matched
            mag_sex_matches_longexp_dict_flat            = flattendict(mag_sex_matches_longexp_dict,ARRAYS=True)
            magerr_sex_matches_longexp_dict_flat         = flattendict(magerr_sex_matches_longexp_dict,ARRAYS=True)
            mag_atlas_matches_longexp_dict_flat          = flattendict(mag_atlas_matches_longexp_dict,ARRAYS=True)
            magerr_atlas_matches_longexp_dict_flat       = flattendict(magerr_atlas_matches_longexp_dict,ARRAYS=True)
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_longexp_dict_flat             = flattendict(zeropoints_all_longexp_dict)
            zeropoints_err_all_longexp_dict_flat         = flattendict(zeropoints_err_all_longexp_dict)
            airmasses_all_longexp_dict_flat              = flattendict(airmasses_all_longexp_dict)
            timeobs_all_longexp_dict_flat                = flattendict(timeobs_all_longexp_dict)
            FWHM_all_longexp_dict_flat                   = flattendict(FWHM_all_longexp_dict)
            donotstackscienceimages_longexp_dict_flat    = flattendict(donotstackscienceimages_longexp_dict)

            ######################## short exposure (120 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_shortexp_dict_flat         = flattendict(stackscienceimages_shortexp_dict)
            zeropoints_shortexp_dict_flat                 = flattendict(zeropoints_shortexp_dict)
            zeropoint_errs_shortexp_dict_flat             = flattendict(zeropoint_errs_shortexp_dict)
            airmasses_shortexp_dict_flat                  = flattendict(airmasses_shortexp_dict)
            timeobs_shortexp_dict_flat                    = flattendict(timeobs_shortexp_dict)
            FWHM_shortexp_dict_flat                       = flattendict(FWHM_shortexp_dict)
            # matched and cleaned
            mag_sex_matches_clean_shortexp_dict_flat      = flattendict(mag_sex_matches_clean_shortexp_dict,ARRAYS=True)
            magerr_sex_matches_clean_shortexp_dict_flat   = flattendict(magerr_sex_matches_clean_shortexp_dict,ARRAYS=True)
            mag_atlas_matches_clean_shortexp_dict_flat    = flattendict(mag_atlas_matches_clean_shortexp_dict,ARRAYS=True)
            magerr_atlas_matches_clean_shortexp_dict_flat = flattendict(magerr_atlas_matches_clean_shortexp_dict,ARRAYS=True)
            # matched
            mag_sex_matches_shortexp_dict_flat            = flattendict(mag_sex_matches_shortexp_dict,ARRAYS=True)
            magerr_sex_matches_shortexp_dict_flat         = flattendict(magerr_sex_matches_shortexp_dict,ARRAYS=True)
            mag_atlas_matches_shortexp_dict_flat          = flattendict(mag_atlas_matches_shortexp_dict,ARRAYS=True)
            magerr_atlas_matches_shortexp_dict_flat       = flattendict(magerr_atlas_matches_shortexp_dict,ARRAYS=True)
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_shortexp_dict_flat             = flattendict(zeropoints_all_shortexp_dict)
            zeropoints_err_all_shortexp_dict_flat         = flattendict(zeropoints_err_all_shortexp_dict)
            airmasses_all_shortexp_dict_flat              = flattendict(airmasses_all_shortexp_dict)
            timeobs_all_shortexp_dict_flat                = flattendict(timeobs_all_shortexp_dict)
            FWHM_all_shortexp_dict_flat                   = flattendict(FWHM_all_shortexp_dict)
            donotstackscienceimages_shortexp_dict_flat    = flattendict(donotstackscienceimages_shortexp_dict)

            if VERBOSE:
                print "Plotting all zero-point data for each passband..."

            # plot zeropoint data for each passband
            plotzp_all(mag_sex_matches_clean_longexp_dict_flat,magerr_sex_matches_clean_longexp_dict_flat,mag_atlas_matches_clean_longexp_dict_flat,magerr_atlas_matches_clean_longexp_dict_flat,mag_sex_matches_longexp_dict_flat,magerr_sex_matches_longexp_dict_flat,mag_atlas_matches_longexp_dict_flat,magerr_atlas_matches_longexp_dict_flat,atlascalibratedplotfolderdir,"LONG",MULTIPROCESSING_CAL,CHECKZP)
            plotzp_all(mag_sex_matches_clean_shortexp_dict_flat,magerr_sex_matches_clean_shortexp_dict_flat,mag_atlas_matches_clean_shortexp_dict_flat,magerr_atlas_matches_clean_shortexp_dict_flat,mag_sex_matches_shortexp_dict_flat,magerr_sex_matches_shortexp_dict_flat,mag_atlas_matches_shortexp_dict_flat,magerr_atlas_matches_shortexp_dict_flat,atlascalibratedplotfolderdir,"SHORT",MULTIPROCESSING_CAL,CHECKZP)

            if VERBOSE:
                print "Plotting zero-point correlations for each passband..."

            # plot zeropoint correlations
            plotzpcorrs(zeropoints_longexp_dict_flat,zeropoint_errs_longexp_dict_flat,zeropoints_all_longexp_dict_flat,zeropoints_err_all_longexp_dict_flat,airmasses_longexp_dict_flat,airmasses_all_longexp_dict_flat,timeobs_longexp_dict_flat,timeobs_all_longexp_dict_flat,FWHM_longexp_dict_flat,FWHM_all_longexp_dict_flat,"LONG",atlascalibratedplotfolderdir,CHECKZP)
            plotzpcorrs(zeropoints_shortexp_dict_flat,zeropoint_errs_shortexp_dict_flat,zeropoints_all_shortexp_dict_flat,zeropoints_err_all_shortexp_dict_flat,airmasses_shortexp_dict_flat,airmasses_all_shortexp_dict_flat,timeobs_shortexp_dict_flat,timeobs_all_shortexp_dict_flat,FWHM_shortexp_dict_flat,FWHM_all_shortexp_dict_flat,"SHORT",atlascalibratedplotfolderdir,CHECKZP)

            # clear variables from memory (running the initial calibration and the second calibration step back-to-back will double their memory and raise a memory error)
            if VERBOSE:
                print "clearing calibration vairables to minimize memory..."

            ######################## long exposure (120 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_longexp_dict = {}             # non-flat dictionaries
            zeropoints_longexp_dict = {}                     # non-flat dictionaries
            zeropoint_errs_longexp_dict = {}                 # non-flat dictionaries
            airmasses_longexp_dict_flat = {}                 # for ZP correlation plots
            airmasses_longexp_dict = {}                      # for ZP correlation plots
            timeobs_longexp_dict_flat = {}                   # for ZP correlation plots
            timeobs_longexp_dict = {}                        # for ZP correlation plots
            FWHM_longexp_dict_flat = {}                      # for ZP correlation plots
            FWHM_longexp_dict = {}                           # for ZP correlation plots
            # matched and cleaned
            mag_sex_matches_clean_longexp_dict_flat = {}     # for ZP plots
            mag_sex_matches_clean_longexp_dict = {}          # for ZP plots
            magerr_sex_matches_clean_longexp_dict_flat = {}  # for ZP plots
            magerr_sex_matches_clean_longexp_dict = {}       # for ZP plots
            mag_ps_matches_clean_longexp_dict_flat = {}      # for ZP plots
            mag_ps_matches_clean_longexp_dict = {}           # for ZP plots
            magerr_ps_matches_clean_longexp_dict_flat = {}   # for ZP plots
            magerr_ps_matches_clean_longexp_dict = {}        # for ZP plots
            # matched
            mag_sex_matches_longexp_dict_flat = {}           # for ZP plots
            mag_sex_matches_longexp_dict = {}                # for ZP plots
            magerr_sex_matches_longexp_dict_flat = {}        # for ZP plots
            magerr_sex_matches_longexp_dict = {}             # for ZP plots
            mag_ps_matches_longexp_dict_flat = {}            # for ZP plots
            mag_ps_matches_longexp_dict = {}                 # for ZP plots
            magerr_ps_matches_longexp_dict_flat = {}         # for ZP plots
            magerr_ps_matches_longexp_dict = {}              # for ZP plots
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_longexp_dict_flat = {}            # for ZP plots
            zeropoints_all_longexp_dict = {}                 # for ZP plots
            zeropoints_err_all_longexp_dict_flat = {}        # for ZP plots
            zeropoints_err_all_longexp_dict = {}             # for ZP plots
            airmasses_all_longexp_dict_flat = {}             # for ZP correlation plots
            airmasses_all_longexp_dict = {}                  # for ZP correlation plots
            timeobs_all_longexp_dict_flat = {}               # for ZP correlation plots
            timeobs_all_longexp_dict = {}                    # for ZP correlation plots
            FWHM_all_longexp_dict_flat = {}                  # for ZP correlation plots
            FWHM_all_longexp_dict = {}                       # for ZP correlation plots
            donotstackscienceimages_longexp_dict = {}        # non-flat dictionary
            ######################## short exposure (5 s) ########################
            # zero slopes and 'good' zero point values
            stackscienceimages_shortexp_dict = {}            # non-flat dictionaries
            zeropoints_shortexp_dict = {}                    # non-flat dictionaries
            zeropoint_errs_shortexp_dict = {}                # non-flat dictionaries
            airmasses_shortexp_dict_flat = {}                # for ZP correlation plots
            airmasses_shortexp_dict = {}                     # for ZP correlation plots
            timeobs_shortexp_dict_flat = {}                  # for ZP correlation plots
            timeobs_shortexp_dict = {}                       # for ZP correlation plots
            FWHM_shortexp_dict_flat = {}                     # for ZP correlation plots
            FWHM_shortexp_dict = {}                          # for ZP correlation plots
            # matched and cleaned
            mag_sex_matches_clean_shortexp_dict_flat = {}    # for ZP plots
            mag_sex_matches_clean_shortexp_dict = {}         # for ZP plots
            magerr_sex_matches_clean_shortexp_dict_flat = {} # for ZP plots
            magerr_sex_matches_clean_shortexp_dict = {}      # for ZP plots
            mag_ps_matches_clean_shortexp_dict_flat = {}     # for ZP plots
            mag_ps_matches_clean_shortexp_dict = {}          # for ZP plots
            magerr_ps_matches_clean_shortexp_dict_flat = {}  # for ZP plots
            magerr_ps_matches_clean_shortexp_dict = {}       # for ZP plots
            # matched
            mag_sex_matches_shortexp_dict_flat = {}          # for ZP plots
            mag_sex_matches_shortexp_dict = {}               # for ZP plots
            magerr_sex_matches_shortexp_dict_flat = {}       # for ZP plots
            magerr_sex_matches_shortexp_dict = {}            # for ZP plots
            mag_ps_matches_shortexp_dict_flat = {}           # for ZP plots
            mag_ps_matches_shortexp_dict = {}                # for ZP plots
            magerr_ps_matches_shortexp_dict_flat = {}        # for ZP plots
            magerr_ps_matches_shortexp_dict = {}             # for ZP plots
            # non-zero slopes and 'bad' zero-point values
            zeropoints_all_shortexp_dict_flat = {}           # for ZP plots
            zeropoints_all_shortexp_dict = {}                # for ZP plots
            zeropoints_err_all_shortexp_dict_flat = {}       # for ZP plots
            zeropoints_err_all_shortexp_dict = {}            # for ZP plots
            airmasses_all_shortexp_dict_flat = {}            # for ZP correlation plots
            airmasses_all_shortexp_dict = {}                 # for ZP correlation plots
            timeobs_all_shortexp_dict_flat = {}              # for ZP correlation plots
            timeobs_all_shortexp_dict = {}                   # for ZP correlation plots
            FWHM_all_shortexp_dict_flat = {}                 # for ZP correlation plots
            FWHM_all_shortexp_dict = {}                      # for ZP correlation plots
            donotstackscienceimages_shortexp_dict = {}       # non-flat dictionary

            print "\nscience images to be stacked (long exp):"
            print stackscienceimages_longexp_dict_flat
            print "\nbad science images (long exposure):"
            print donotstackscienceimages_longexp_dict_flat
            print "\nscience images to be stacked (short exp):"
            print stackscienceimages_shortexp_dict_flat
            print "\nbad science images (short exposure):"
            print donotstackscienceimages_shortexp_dict_flat
            #print "\nlogstring:"
            #print logstring
        else:
            if VERBOSE:
                print "Empty science files:", sciencefiles
                print "Skipping calibration check."

        if LOGOUTPUT:
            # write calibration log
            if VERBOSE:
                print "writing zeropoint log"
            f=open(logdir+"zeropoint.log","w+")
            for item in logstring:
                f.write("%s\n" % item)
            f.close()


#########################################################################################################
#                                                                                                       #
#                                             STACKING                                                  #
#                                                                                                       #
#########################################################################################################

if ATLAS and STACK:
    CHECKZP=True

    if VERBOSE:
        print "\n-------------------- proceeding with stacking calibrated images ------------------"

    # clear catalog directories
    if pathexists(stackedcaldir):
        os.system("rm -rf "+stackedcaldir)
        os.system("mkdir "+stackedcaldir)
        if VERBOSE:
            print "cleared directory "+stackedcaldir
    else:
        os.system("mkdir "+stackedcaldir)
        if VERBOSE:
            print "created directory "+stackedcaldir

    stackedcalframesdir = stackedcaldir+"frames/"

    os.system("mkdir "+stackedcalframesdir)

    # plots
    stackedcalplotsdir  = stackedcaldir+"plots/"
    if pathexists(stackedcalplotsdir):
        os.system("rm -rf "+stackedcalplotsdir)
        if VERBOSE:
            print "cleared directory "+stackedcalplotsdir
    os.system("mkdir "+stackedcalplotsdir)
    if VERBOSE:
        print "created directory "+stackedcalplotsdir

    #################################### long exposure (120s) ####################################

    stackscienceimages_longexp_dict = freaddict_pickle("stackscienceimages_all_longexp_dict_pickle.txt")
    
    if VERBOSE:
        print "stackscienceimages_longexp_dict:"
        print stackscienceimages_longexp_dict

    # create sub-folders for each pointing
    for object in stackscienceimages_longexp_dict.keys():
        # images
        os.system("mkdir "+stackedcalframesdir+object)
        if VERBOSE:
            print "created directory "+stackedcalframesdir+object
        # plots
        os.system("mkdir "+stackedcalplotsdir+object+"/")
        if VERBOSE:
            print " created directory "+stackedcalplotsdir+object+"/"
    
    if VERBOSE:
        print "\n----------------------------- stacking long exposures ----------------------------"

    # main median-stacking function
    fcalstack_main(stackscienceimages_longexp_dict,stackedcalframesdir,"LONG")

    #################################### short exposure (5s) #####################################

    stackscienceimages_shortexp_dict = freaddict_pickle("stackscienceimages_all_shortexp_dict_pickle.txt")
    
    if VERBOSE:
       print "stackscienceimages_shortexp_dict:"
       print stackscienceimages_shortexp_dict

    # create sub-folders for each pointing
    for object in stackscienceimages_shortexp_dict.keys():
        if not pathexists(stackedcalframesdir+object):
            # images
            os.system("mkdir "+stackedcalframesdir+object)
            if VERBOSE:
                print "created directory "+stackedcalframesdir+object

    if VERBOSE:
        print "\n----------------------------- stacking short exposures ---------------------------"

    # main median-stacking function
    fcalstack_main(stackscienceimages_shortexp_dict,stackedcalframesdir,"SHORT")

else:
    if VERBOSE:
        print "\n----------------------- skipping stacking calibrated images ----------------------"

#########################################################################################################
#                                                                                                       #
#                                      FINAL CALIBRATION CHECK                                          #
#                                                                                                       #
#########################################################################################################

# zero-point of median-stacked frames should be consistent with that of individual frames

if ATLAS and FINALCALCHECK:
    FINALCHECKZP = True

    if VERBOSE:
        print "\n---------------- proceeding with calibration check on stacked images -------------"

    atlastrimcatf     = atlasdir+atlastrimcat
    hdul_atlastrim    = fits.open(atlastrimcatf)
    atlasdata_trim    = hdul_atlastrim[1].data

    stackedcalframesdir = stackedcaldir+"frames/"
    stackedcalplotsdir  = stackedcaldir+"plots/"

    # SExtractor folder
    if pathexists(sexfinalcheckzpdir):
        os.system("rm -rf "+sexfinalcheckzpdir)
        if VERBOSE:
            print "cleared directory "+sexfinalcheckzpdir
    os.system("mkdir "+sexfinalcheckzpdir)
    if VERBOSE:
        print "created directory "+sexfinalcheckzpdir

    # objects
    os.system("mkdir "+sexobjfinalcheckzpdir)
    if VERBOSE:
        print "created directory "+sexobjfinalcheckzpdir
    # catalogues
    os.system("mkdir "+sexcatfinalcheckzpdir)
    if VERBOSE:
        print "created directory "+sexcatfinalcheckzpdir
    # background
    os.system("mkdir "+sexbgfinalcheckzpdir)
    if VERBOSE:
        print "created directory "+sexbgfinalcheckzpdir

    # plots
    stackedcalplotsdir  = stackedcaldir+"plots/"
    if not pathexists(stackedcalplotsdir):
        os.system("mkdir "+stackedcalplotsdir)
        if VERBOSE:
            print "created directory "+stackedcalplotsdir

    #################################### long exposure (120s) ####################################

    stackscienceimages_longexp_dict = freaddict_pickle("stackscienceimages_all_longexp_dict_pickle.txt")

    # Source Extractor
    for object in stackscienceimages_longexp_dict.keys():
        # objects
        os.system("mkdir "+sexobjfinalcheckzpdir+object+"/")
        # catalogues
        os.system("mkdir "+sexcatfinalcheckzpdir+object+"/")
        # background
        os.system("mkdir "+sexbgfinalcheckzpdir+object+"/")
        # plots
        #os.system("mkdir "+stackedcalplotsdir+object+"/")

    objects = stackscienceimages_longexp_dict.keys()

    for object in objects:
        if  VERBOSE:
            print "processing long exposures of pointing "+object+"..."

        zeropoints_longexp_r   = []
        zeropoints_longexp_i   = []
        zeropoints_longexp_z   = []

        stackedcalframesobjdir = stackedcalframesdir+object+"/"
        sciencefiles           = os.listdir(stackedcalframesobjdir)
        for file in sciencefiles:
            _=mainzeropointcal_atlas(atlasdata_trim,object,stackedcalframesobjdir,file,logbool=LOGOUTPUT,FINALCHECKZP=FINALCHECKZP)

    #################################### short exposure (5s) #####################################

    stackscienceimages_shortexp_dict = freaddict_pickle("stackscienceimages_all_shortexp_dict_pickle.txt")

    objects = stackscienceimages_shortexp_dict.keys()

    for object in objects:
        if  VERBOSE:
            print "processing short exposures of pointing "+object+"..."

        zeropoints_shortexp_r  = []

        stackedcalframesobjdir = stackedcalframesdir+object+"/"
        sciencefiles = os.listdir(stackedcalframesobjdir)
        for file in sciencefiles:
            _=mainzeropointcal_atlas(atlasdata_trim,object,stackedcalframesobjdir,file,logbool=LOGOUTPUT,FINALCHECKZP=FINALCHECKZP)
    
else:
    if VERBOSE:
        print "\n------------------- skipping calibration check on stacked images ----------------"

#########################################################################################################
#                                                                                                       #
#                                      FINAL SOURCE EXTRACTOR RUN                                       #
#                                                                                                       #
#########################################################################################################

if ATLAS and FINALSEXRUN:
    FINALCHECKZP = True

    if VERBOSE:
        print "\n--------------------- creating final Source Extractor catalogues -----------------"

    stackedcalframesdir = stackedcaldir+"frames/"

    # SExtractor folder
    if pathexists(sexstackedcaldir):
        os.system("rm -rf "+sexstackedcaldir)
        if VERBOSE:
            print "cleared directory "+sexstackedcaldir
    os.system("mkdir "+sexstackedcaldir)
    if VERBOSE:
        print "created directory "+sexstackedcaldir

    # objects
    os.system("mkdir "+sexobjstackedcaldir)
    if VERBOSE:
        print "created directory "+sexobjstackedcaldir
    # catalogues
    os.system("mkdir "+sexcatstackedcaldir)
    if VERBOSE:
        print "created directory "+sexcatstackedcaldir
    # background
    os.system("mkdir "+sexbgstackedcaldir)
    if VERBOSE:
        print "created directory "+sexbgstackedcaldir

    for object in os.listdir(stackedcalframesdir):
        # objects
        os.system("mkdir "+sexobjstackedcaldir+object+"/")
        # catalogues
        os.system("mkdir "+sexcatstackedcaldir+object+"/")
        # background
        os.system("mkdir "+sexbgstackedcaldir+object+"/")

    for object in os.listdir(stackedcalframesdir):
        if VERBOSE:
            print "processing pointing "+object+"..."
        for f in os.listdir(stackedcalframesdir+object+"/"):
            fpath = stackedcalframesdir+object+"/"
            #sexcall(f,object,fpath,sexobjstackedcaldir,sexcatstackedcaldir,sexbgstackedcaldir,sexbool=SEXTRACTOR,FINALCHECKZP=FINALCHECKZP)
            # going to brute-force this to make sure we add zero-point to call for proper calibration this time (no longer inspecting zero-point)
            fname         = f.split(".fts")[0]
            header        = fits.getheader(fpath+f)
            MAG_ZEROPOINT = header["ZP"]
            # Construct source extractor calls
            objsexcall = "sex -MAG_ZEROPOINT "+str(MAG_ZEROPOINT)+" -CATALOG_TYPE FITS_1.0 -PARAMETERS_NAME default.param -CATALOG_NAME "+sexcatstackedcaldir+object+"/"+fname+".fts"+" -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME "+sexobjstackedcaldir+object+"/"+fname+"_objects.fts "+fpath+f
            baksexcall = "sex -MAG_ZEROPOINT "+str(MAG_ZEROPOINT)+" -CATALOG_TYPE FITS_1.0 -PARAMETERS_NAME default.param -CATALOG_NAME "+sexcatstackedcaldir+object+"/"+fname+".fts"+" -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME "+sexbgstackedcaldir+object+"/"+fname+"_background.fts "+fpath+f
            os.system(objsexcall)
            os.system(baksexcall)
            
else:
    if VERBOSE:
        print "\n--------------------- skipping final Source Extractor catalogue -----------------"


#########################################################################################################
#                                                                                                       #
#                                      MATCH SOURCE CATALOGUES                                          #
#                                                                                                       #
#########################################################################################################

ATLASDITcatf = mergedcatsdir+ATLASDITcat # merged catalogue between ATLAS and DIT
if ATLAS and MATCHCATS:
    if VERBOSE:
        print "\n------------------------ Matching Source Extractor catalogues -------------------"

    # successively merge DIT catalogues with merged ATLAS catalogue

    objects             = ["W3_1","W3_11","W3_12","W3_13","W3_21","W3_22","W3_23","W3_31","W3_32","W3_33"]

    #################################### long exposure (120s) ####################################

    if VERBOSE:
        print "Merging long (120 s) exposures..."

    stackscienceimages_longexp_dict = freaddict_pickle("stackscienceimages_all_longexp_dict_pickle.txt")

    for object in objects:
        stackedcalcatobjdir = sexcatstackedcaldir+object+"/"
        if pathexists(stackedcalcatobjdir):
            if VERBOSE:
                print ">>> working on directory: "+stackedcalcatobjdir
            files = os.listdir(stackedcalcatobjdir)
            for file in files:
                exp = file.split(".cat")[0].split("_")[-1]
                if exp=="longexp":
                    DITcatf = stackedcalcatobjdir+file
                    stiltscall = "stilts tskymatch2 in1={} in2={} out={} ra1=ra dec1=dec ra2=RA dec2=Dec error=1 find=best".format(ATLASDITcatf,DITcatf,ATLASDITcatf)
                    print stiltscall
        else:
            ">>> directory empty: "+stackedcalcatobjdir

    #################################### short exposure (5s) #####################################

    stackscienceimages_shortexp_dict = freaddict_pickle("stackscienceimages_all_shortexp_dict_pickle.txt")

    if VERBOSE:
        print "Merging short (5 s) exposures..."

    for object in objects:
        stackedcalcatobjdir = sexcatstackedcaldir+object+"/"
        if pathexists(stackedcalcatobjdir):
            if VERBOSE:
                print ">>> working on directory: "+stackedcalcatobjdir
            files = os.listdir(stackedcalcatobjdir)
            for file in files:
                exp = file.split(".cat")[0].split("_")[-1]
                if exp=="shortexp":
                    DITcatf = stackedcalcatobjdir+file
                    stiltscall = "stilts tskymatch2 in1={} in2={} out={} ra1=ra dec1=dec ra2=RA dec2=Dec error=1 find=best".format(ATLASDITcatf,DITcatf,ATLASDITcatf)
                    print stiltscall
        else:
            ">>> directory empty: "+stackedcalcatobjdir


    #atlascatf  = atlasdir+atlascat
    #stiltscall = "stilts tskymatch2 in1={} in2={} out={} ra1=ra dec1=dec ra2=RA dec2=Dec error=1 find=best".format(atlascatf,DITcatf,atlasDITcatf)


else:
    if VERBOSE:
        print "\n-------------------- skipping Source Extractor catalogue matching ---------------"



