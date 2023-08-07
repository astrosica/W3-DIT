#!/usr/local/bin/python

import os
import yaml
import pylab
import numpy as np 
from astropy.io import fits

import matplotlib.pyplot as plt

#########################################################################################################
#                                                                                                       #
#                                 IMPORT CONFIGURATION FILE                                             #
#                                                                                                       #
#########################################################################################################

with open("config.yml","r") as fconfig:
    config_data = yaml.safe_load(fconfig)

# directories
parentdir       = config_data["parentdir"]
mergedcatsdir   = config_data["mergedcatsdir"]
ATLAS_DIT_cat   = config_data["ATLASDITcat"]
plotdir         = parentdir+"plots/compare-DIT-ATLAS/"

#########################################################################################################
#                                                                                                       #
#                                          DEFINE FUNCTIONS                                             #
#                                                                                                       #
#########################################################################################################

def pathexists(path):
    '''
    path: directory to check existence for
    Returns True if path exists, False if it does not
    '''
    if os.path.exists(path):
        return True
    else:
        return False

def fmagdiff(mag_DIT,magerr_DIT,mag_ATLAS,magerr_ATLAS):
    '''
    Calculates (DIT mag - ATLAS mag) and associated uncertainties.
    '''
    
    magdiff        = np.array(mag_ATLAS) - np.array(mag_DIT)
    magdiff_error  = np.sqrt((magerr_ATLAS)**2. + (magerr_DIT)**2.)
    
    return magdiff,magdiff_error

def fplotfigs(DIT_mag,DIT_magerr,ATLAS_mag,ATLAS_magerr,passband,fpath):
	######### DIT VS ATLAS #########
	# scatterplot
	fig = plt.figure(1,figsize=(11,8.5))
	plt.errorbar(DIT_mag,ATLAS_mag,xerr=DIT_magerr,yerr=ATLAS_magerr,fmt="o",markersize=8,color="red")
	plt.xlabel("ATLAS {} [mag]".format(passband))
	plt.ylabel("DIT {} [mag]".format(passband,passband))
	plt.savefig(fpath+"_corr.png",bbox_inches="tight")
	plt.close()
	# 2d histogram
	fig = plt.figure(1,figsize=(11,8.5))
	plt.hist2d(DIT_mag,ATLAS_mag,bins=(100,100),cmin=1,cmap=plt.cm.plasma)
	plt.xlabel("ATLAS {} [mag]".format(passband))
	plt.ylabel("DIT {} [mag]".format(passband,passband))
	print "saving figure "+fpath+"_corr_2dhist.png"
	plt.savefig(fpath+"_corr_2dhist.png",bbox_inches="tight")
	plt.close()

	magdiff,magdiff_err = fmagdiff(DIT_mag,DIT_magerr,ATLAS_mag,ATLAS_magerr)
	##### (DIT-ATLAS) vs ATLAS #####
	# scatterplot
	fig = plt.figure(1,figsize=(11,8.5))
	plt.errorbar(ATLAS_mag,magdiff,xerr=ATLAS_magerr,yerr=magdiff_err,fmt="o",markersize=8,color="red")
	plt.xlabel("ATLAS {} [mag]".format(passband))
	plt.ylabel("(DIT {} - ATLAS {}) [mag]".format(passband,passband))
	plt.savefig(fpath+"_zp.png",bbox_inches="tight")
	plt.close()
	# 2d histogram
	fig = plt.figure(1,figsize=(11,8.5))
	plt.hist2d(ATLAS_mag,magdiff,bins=(100,70),cmin=1,cmap=plt.cm.plasma)
	plt.xlabel("ATLAS {} [mag]".format(passband))
	plt.ylabel("(DIT {} - ATLAS {}) [mag]".format(passband,passband))
	print "saving figure "+fpath+"_zp_2dhist.png"
	plt.savefig(fpath+"_zp_2dhist.png",bbox_inches="tight")
	plt.close()

def fplotmaghist(mags,passband,fpath):
	'''
	'''
	#N      = len(np.arange(poslower,posupper,deltapos))
	#params = dict(bins=N,range=(poslower,posupper))
	fig = plt.figure(1,figsize=(11,8.5))
	ax  = fig.add_subplot(111)
	pylab.hist(mags)#,**params)
	#plt.axvline(posmatch,linestyle="dashed",color="black",label="position match")
	plt.xlabel("DIT {} [mags]".format(passband))
	plt.ylabel("N")
	print "saving figure "+fpath+"_DITmags_noATLAS_hist.png"
	plt.savefig(fpath+"_DITmags_noATLAS_hist.png",bbox_inches="tight")
	plt.close()


#########################################################################################################
#                                                                                                       #
#                                    READ CATALOGUE DATA                                                #
#                                                                                                       #
#########################################################################################################

ATLAS_DIT_catf = mergedcatsdir+ATLAS_DIT_cat

print "Reading merged catalogue..."

hdul = fits.open(ATLAS_DIT_catf)
data = hdul[1].data

########### ATLAS ###########

print "Extrtacting ATLAS catalogue..."

RA_ATLAS   = data["RA"]
DEC_ATLAS  = data["DEC"]
r_ATLAS    = data["r"]
rerr_ATLAS = data["dr"]
i_ATLAS    = data["i"]
ierr_ATLAS = data["di"]
z_ATLAS    = data["z"]
zerr_ATLAS = data["dz"]

############ DIT ############

# LONG (120s) EXPOSURES

print "Extracting DIT long exposure (120s) catalogues..."

# pointing 1 (_1): r-band
RA_DIT_1   = data["X_WORLD_1"]
DEC_DIT_1  = data["Y_WORLD_1"]
r_DIT_1    = data["MAG_AUTO_1"]
rerr_DIT_1 = data["MAGERR_AUTO_1"]
# pointing 2 (_2): i-band
RA_DIT_2   = data["X_WORLD_2"]
DEC_DIT_2  = data["Y_WORLD_2"]
i_DIT_2    = data["MAG_AUTO_2"]
ierr_DIT_2 = data["MAGERR_AUTO_2"]
# pointing 3 (_1a): z-band
RA_DIT_3   = data["X_WORLD_1a"]
DEC_DIT_3  = data["Y_WORLD_1a"]
z_DIT_3    = data["MAG_AUTO_1a"]
zerr_DIT_3 = data["MAGERR_AUTO_1a"]
# pointing 4 (_2a): i-band
RA_DIT_4   = data["X_WORLD_2a"]
DEC_DIT_4  = data["Y_WORLD_2a"]
i_DIT_4    = data["MAG_AUTO_2a"]
ierr_DIT_4 = data["MAGERR_AUTO_2a"]
# pointing 5 (_1b): r-band
RA_DIT_5   = data["X_WORLD_1b"]
DEC_DIT_5  = data["Y_WORLD_1b"]
r_DIT_5    = data["MAG_AUTO_1b"]
rerr_DIT_5 = data["MAGERR_AUTO_1b"]
# pointing 6 (_2b): z-band
RA_DIT_6   = data["X_WORLD_2b"]
DEC_DIT_6  = data["Y_WORLD_2b"]
z_DIT_6    = data["MAG_AUTO_2b"]
zerr_DIT_6 = data["MAGERR_AUTO_2b"]
# pointing 7 (_1c): r-band
RA_DIT_7   = data["X_WORLD_1c"]
DEC_DIT_7  = data["Y_WORLD_1c"]
r_DIT_7    = data["MAG_AUTO_1c"]
rerr_DIT_7 = data["MAGERR_AUTO_1c"]
# pointing 8 (_2c): i-band
RA_DIT_8   = data["X_WORLD_2c"]
DEC_DIT_8  = data["Y_WORLD_2c"]
i_DIT_8    = data["MAG_AUTO_2c"]
ierr_DIT_8 = data["MAGERR_AUTO_2c"]
# pointing 9 (_1d): z-band
RA_DIT_9   = data["X_WORLD_1d"]
DEC_DIT_9  = data["Y_WORLD_1d"]
z_DIT_9    = data["MAG_AUTO_1d"]
zerr_DIT_9 = data["MAGERR_AUTO_1d"]
# pointing 10 (_2d): z-band
RA_DIT_10   = data["X_WORLD_2d"]
DEC_DIT_10  = data["Y_WORLD_2d"]
z_DIT_10    = data["MAG_AUTO_2d"]
zerr_DIT_10 = data["MAGERR_AUTO_2d"]
# pointing 11 (_1e): i-band
RA_DIT_11   = data["X_WORLD_1e"]
DEC_DIT_11  = data["Y_WORLD_1e"]
i_DIT_11    = data["MAG_AUTO_1e"]
ierr_DIT_11 = data["MAGERR_AUTO_1e"]
# pointing 12 (_2e): r-band
RA_DIT_12   = data["X_WORLD_2e"]
DEC_DIT_12  = data["Y_WORLD_2e"]
r_DIT_12    = data["MAG_AUTO_2e"]
rerr_DIT_12 = data["MAGERR_AUTO_2e"]
# pointing 13 (_1f): i-band
RA_DIT_13   = data["X_WORLD_1f"]
DEC_DIT_13  = data["Y_WORLD_1f"]
i_DIT_13    = data["MAG_AUTO_1f"]
ierr_DIT_13 = data["MAGERR_AUTO_1f"]
# pointing 14 (_2f): r-band
RA_DIT_14   = data["X_WORLD_2f"]
DEC_DIT_14  = data["Y_WORLD_2f"]
r_DIT_14    = data["MAG_AUTO_2f"]
rerr_DIT_14 = data["MAGERR_AUTO_2f"]
# pointing 15 (_1g): z-band
RA_DIT_15   = data["X_WORLD_1g"]
DEC_DIT_15  = data["Y_WORLD_1g"]
z_DIT_15    = data["MAG_AUTO_1g"]
zerr_DIT_15 = data["MAGERR_AUTO_1g"]
# pointing 16 (_2g): i-band
RA_DIT_16   = data["X_WORLD_2g"]
DEC_DIT_16  = data["Y_WORLD_2g"]
i_DIT_16    = data["MAG_AUTO_2g"]
ierr_DIT_16 = data["MAGERR_AUTO_2g"]
# pointing 17 (_1h): r-band
RA_DIT_17   = data["X_WORLD_1h"]
DEC_DIT_17  = data["Y_WORLD_1h"]
r_DIT_17    = data["MAG_AUTO_1h"]
rerr_DIT_17 = data["MAGERR_AUTO_1h"]
# pointing 18 (_2h): z-band
RA_DIT_18   = data["X_WORLD_2h"]
DEC_DIT_18  = data["Y_WORLD_2h"]
z_DIT_18    = data["MAG_AUTO_2h"]
zerr_DIT_18 = data["MAGERR_AUTO_2h"]
# pointing 19 (_1i): r-band
RA_DIT_19   = data["X_WORLD_1i"]
DEC_DIT_19  = data["Y_WORLD_1i"]
r_DIT_19    = data["MAG_AUTO_1i"]
rerr_DIT_19 = data["MAGERR_AUTO_1i"]
# pointing 20 (_2i): z-band
RA_DIT_20   = data["X_WORLD_2i"]
DEC_DIT_20  = data["Y_WORLD_2i"]
z_DIT_20    = data["MAG_AUTO_2i"]
zerr_DIT_20 = data["MAGERR_AUTO_2i"]

# SHORT (5s) EXPOSURES

print "Extracting DIT short exposure (5s) catalogues..."

# pointing 21 (_1j): r-band
RA_DIT_21   = data["X_WORLD_1j"]
DEC_DIT_21  = data["Y_WORLD_1j"]
r_DIT_21    = data["MAG_AUTO_1j"]
rerr_DIT_21 = data["MAGERR_AUTO_1j"]
# pointing 22 (_2j): r-band
RA_DIT_22   = data["X_WORLD_2j"]
DEC_DIT_22  = data["Y_WORLD_2j"]
r_DIT_22    = data["MAG_AUTO_2j"]
rerr_DIT_22 = data["MAGERR_AUTO_2j"]
# pointing 23 (_1k): r-band
RA_DIT_23   = data["X_WORLD_1k"]
DEC_DIT_23  = data["Y_WORLD_1k"]
r_DIT_23    = data["MAG_AUTO_1k"]
rerr_DIT_23 = data["MAGERR_AUTO_1k"]
# pointing 24 (_2k): r-band
RA_DIT_24   = data["X_WORLD_2k"]
DEC_DIT_24  = data["Y_WORLD_2k"]
r_DIT_24    = data["MAG_AUTO_2k"]
rerr_DIT_24 = data["MAGERR_AUTO_2k"]

#########################################################################################################
#                                                                                                       #
#                                    COMPARE DIT TO ATLAS                                               #
#                                                                                                       #
#########################################################################################################

if pathexists(plotdir):
	os.system("rm -rf "+plotdir)
	print "cleared directory "+plotdir
os.system("mkdir "+plotdir)
print "created directory "+plotdir

pointings = ["W3_1","W3_11","W3_12","W3_13","W3_31","W3_32","W3_33"]
for pointing in pointings:
	os.system("mkdir "+plotdir+pointing+"/")
	print "created directory "+plotdir+pointing+"/"

# LONG (120s) EXPOSURES
print "Plotting long exposure (120s) figures..."
exp = "longexp"

# pointing 1 (_1): r-band

passband            = "r"
pointing            = "W3_1"

DIT_mag             = r_DIT_1
DIT_magerr          = rerr_DIT_1
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 2 (_2): i-band

passband            = "i"
pointing            = "W3_1"

DIT_mag             = i_DIT_2
DIT_magerr          = ierr_DIT_2
ATLAS_mag           = i_ATLAS
ATLAS_magerr        = ierr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 3 (_1a): z-band

passband            = "z"
pointing            = "W3_1"

DIT_mag             = z_DIT_3
DIT_magerr          = zerr_DIT_3
ATLAS_mag           = z_ATLAS
ATLAS_magerr        = zerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 4 (_2a): i-band

passband            = "i"
pointing            = "W3_11"

DIT_mag             = i_DIT_4
DIT_magerr          = ierr_DIT_4
ATLAS_mag           = i_ATLAS
ATLAS_magerr        = ierr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 5 (_1b): r-band

passband            = "r"
pointing            = "W3_11"

DIT_mag             = r_DIT_5
DIT_magerr          = rerr_DIT_5
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 6 (_2b): z-band

passband            = "z"
pointing            = "W3_11"

DIT_mag             = z_DIT_6
DIT_magerr          = zerr_DIT_6
ATLAS_mag           = z_ATLAS
ATLAS_magerr        = zerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 7 (_1c): r-band

passband            = "r"
pointing            = "W3_12"

DIT_mag             = r_DIT_7
DIT_magerr          = rerr_DIT_7
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 8 (_2c): i-band

passband            = "i"
pointing            = "W3_12"

DIT_mag             = i_DIT_8
DIT_magerr          = ierr_DIT_8
ATLAS_mag           = i_ATLAS
ATLAS_magerr        = ierr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 9 (_1d): z-band

passband            = "z"
pointing            = "W3_12"

DIT_mag             = z_DIT_9
DIT_magerr          = zerr_DIT_9
ATLAS_mag           = z_ATLAS
ATLAS_magerr        = zerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 10 (_2d): z-band

passband            = "z"
pointing            = "W3_13"

DIT_mag             = z_DIT_10
DIT_magerr          = zerr_DIT_10
ATLAS_mag           = z_ATLAS
ATLAS_magerr        = zerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 11 (_1e): i-band

passband            = "i"
pointing            = "W3_13"

DIT_mag             = i_DIT_11
DIT_magerr          = ierr_DIT_11
ATLAS_mag           = i_ATLAS
ATLAS_magerr        = ierr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 12 (_2e): r-band

passband            = "r"
pointing            = "W3_13"

DIT_mag             = r_DIT_12
DIT_magerr          = rerr_DIT_12
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 13 (_1f): i-band

passband            = "i"
pointing            = "W3_31"

DIT_mag             = i_DIT_13
DIT_magerr          = ierr_DIT_13
ATLAS_mag           = i_ATLAS
ATLAS_magerr        = ierr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 14 (_2f): r-band

passband            = "r"
pointing            = "W3_31"

DIT_mag             = r_DIT_14
DIT_magerr          = rerr_DIT_14
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 15 (_1g): z-band

passband            = "z"
pointing            = "W3_31"

DIT_mag             = z_DIT_15
DIT_magerr          = zerr_DIT_15
ATLAS_mag           = z_ATLAS
ATLAS_magerr        = zerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 16 (_2g): i-band

passband            = "i"
pointing            = "W3_32"

DIT_mag             = i_DIT_16
DIT_magerr          = ierr_DIT_16
ATLAS_mag           = i_ATLAS
ATLAS_magerr        = ierr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 17 (_1h): r-band

passband            = "r"
pointing            = "W3_32"

DIT_mag             = r_DIT_17
DIT_magerr          = rerr_DIT_17
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 18 (_2h): z-band

passband            = "z"
pointing            = "W3_32"

DIT_mag             = z_DIT_18
DIT_magerr          = zerr_DIT_18
ATLAS_mag           = z_ATLAS
ATLAS_magerr        = zerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 19 (_1i): r-band

passband            = "r"
pointing            = "W3_33"

DIT_mag             = r_DIT_19
DIT_magerr          = rerr_DIT_19
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 20 (_2i): z-band

passband            = "z"
pointing            = "W3_33"

DIT_mag             = z_DIT_20
DIT_magerr          = zerr_DIT_20
ATLAS_mag           = z_ATLAS
ATLAS_magerr        = zerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)




# SHORT (5s) EXPOSURES
print "Plotting long exposure (120s) figures..."
exp      = "shortexp"
passband = "r"

# pointing 21 (_1j): r-band

pointing            = "W3_1"

DIT_mag             = r_DIT_21
DIT_magerr          = rerr_DIT_21
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 22 (_2j): r-band

pointing            = "W3_11"

DIT_mag             = r_DIT_22
DIT_magerr          = rerr_DIT_22
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 23 (_1k): r-band

pointing            = "W3_13"

DIT_mag             = r_DIT_23
DIT_magerr          = rerr_DIT_23
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)

# pointing 24 (_2k): r-band

pointing            = "W3_32"

DIT_mag             = r_DIT_24
DIT_magerr          = rerr_DIT_24
ATLAS_mag           = r_ATLAS
ATLAS_magerr        = rerr_ATLAS

# remove nans
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==False) )[0]
DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean = DIT_mag[ii],DIT_magerr[ii],ATLAS_mag[ii],ATLAS_magerr[ii]

# plot
pointing_plotdir = plotdir+pointing+"/"
fname            = "{}_{}_{}".format(pointing,passband,exp)
fplotfigs(DIT_mag_clean,DIT_magerr_clean,ATLAS_mag_clean,ATLAS_magerr_clean,passband,pointing_plotdir+fname)

# DIT sources with no ATLAS counterpart
ii = np.where( (np.isnan(DIT_mag)==False) & (np.isnan(ATLAS_mag)==True) )[0]
DIT_mag_noATLAS     = DIT_mag[ii]
DIT_magerr_noATLAS  = DIT_magerr[ii]

if len(DIT_mag_noATLAS)>0:
	fplotmaghist(DIT_mag_noATLAS,passband,pointing_plotdir+fname)




