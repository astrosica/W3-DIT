#!/usr/local/bin/python

import pylab
import numpy as np 
from scipy import spatial
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import matplotlib.pyplot as plt

#########################################################################################################
#                                                                                                       #
#                                 IMPORT CONFIGURATION FILE                                             #
#                                                                                                       #
#########################################################################################################

atlasdir      = "/mnt/raid-project/hp/campbell/ATLAS/"
atlascatf     = "ATLASw3_pgmartin.fits"
panstarrsdir  = "/mnt/raid-project/hp/campbell/panstarrs/"
panstarrscatf = "ps_box.fits"
panstarrsphot = "PSF"
posmatch      = 1.0 # arcsecond
plotdir       = "/mnt/raid-project/hp/campbell/DIT/Python/Pipeline/plots/atlas-panstarrs/"

#########################################################################################################
#                                                                                                       #
#                                        FUNCTIONS                                                      #
#                                                                                                       #
#########################################################################################################

def pointingoffsethist(file,ra_atlas,dec_atlas,ra_ps,dec_ps,indices,dist,plotdir,radius=posmatch,deltapos=1.0,poslower=0.0,posupper=10.0):
    '''
    Plots the pointing offset between ATLAS and Source Extractor catalog matches from KDTree.
    '''
    
    ra_atlas_matches   = ra_atlas[indices]
    dec_atlas_matches  = dec_atlas[indices]
    
    deltara           = abs(np.array(ra_atlas_matches-ra_ps))
    deltadec          = abs(np.array(dec_atlas_matches-dec_ps))
    deltapos_deg      = np.sqrt((deltara)**2. + (deltadec)**2.)
    deltapos_arcsec   = deltapos_deg*3600.
    dist_arcsec       = dist*3600.
    
    newfname = file+"_posresidualhist.pdf"
    
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
    plt.savefig(plotdir+newfname,bbox_inches="tight")
    plt.close(fig)

#########################################################################################################
#                                                                                                       #
#                                    Pan-STARRS CATALOGUE                                               #
#                                                                                                       #
#########################################################################################################

print "Reading in Pan-STARRS catalogue..."
hdul_panstarrs       = fits.open(panstarrsdir+panstarrscatf)
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


#########################################################################################################
#                                                                                                       #
#                                       ATLAS CATALOGUE                                                 #
#                                                                                                       #
#########################################################################################################

print "Reading in ATLAS catalogue..."
hdul_atlas       = fits.open(atlasdir+atlascatf)
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

#########################################################################################################
#                                                                                                       #
#                                      MATCH CATALOGUES                                                 #
#                                                                                                       #
#########################################################################################################

print "Matching ATLAS and Pan-STARRS..."

# usage        : match_coordinates_sky(matchcoord, catalogcoord, nthneighbor=1, storekdtree='kdtree_sky')
#
# Inputs
# matchcoord   : The coordinate(s) to match to the catalog.
# catalogcoord : The base catalog in which to search for matches. Typically this will be a coordinate object that is an array.
# nthneighbor  : Which closest neighbor to search for. Typically 1 is desired here, as that is correct for matching one set of coordinates to another. 
# 
# Outputs
# idx          : Indices into catalogcoord to get the matched points for each matchcoord. Shape matches matchcoord.
# sep2d        : The on-sky separation between the closest match for each matchcoord and the matchcoord. Shape matches matchcoord.
# dist3d       : The 3D distance between the closest match for each matchcoord and the matchcoord. Shape matches matchcoord.

# search ATLAS to match to Pan-STARRS (because Pan-STARRS is a much bigger catalogue)
PS_coord                   = SkyCoord(ra=PS_raMean*u.degree, dec=PS_decMean*u.degree)
ATLAS_coord                = SkyCoord(ra=ATLAS_RA*u.degree,  dec=ATLAS_DEC*u.degree)
idx_closest,d2d_closest,_  = match_coordinates_sky(PS_coord,ATLAS_coord,nthneighbor=1)
d2d_arcsec_closest         = d2d_closest.arcsecond

N = len(np.arange(0.,10.,0.1))
params = dict(bins=N,range=(0.,10.))

fig = plt.figure(1,figsize=(11,8.5))
ax = fig.add_subplot(111)
pylab.hist(d2d_arcsec_closest,**params)
plt.axvline(1.0,linestyle="dashed",color="black",label="position match")
plt.xlabel("Residual Pointing Offset (arcseconds)")
plt.ylabel("N")
plt.legend(loc="upper right")
plt.savefig(plotdir+"panstarrs_vs_atlas_posresidualhist.pdf",bbox_inches="tight")
plt.close(fig)

# ATLAS: closest to Pan-STARRS
ATLAS_RA_closest       = ATLAS_RA[idx_closest]
ATLAS_DEC_closest      = ATLAS_DEC[idx_closest]
ATLAS_g_closest        = ATLAS_g[idx_closest]
ATLAS_dg_closest       = ATLAS_dg[idx_closest]
ATLAS_r_closest        = ATLAS_r[idx_closest]
ATLAS_dr_closest       = ATLAS_dr[idx_closest]
ATLAS_i_closest        = ATLAS_i[idx_closest]
ATLAS_di_closest       = ATLAS_di[idx_closest]
ATLAS_z_closest        = ATLAS_z[idx_closest]
ATLAS_dz_closest       = ATLAS_dz[idx_closest]

#########################################################################################################
#                                                                                                       #
#                                      CLEANN MATCHING                                                  #
#                                                                                                       #
#########################################################################################################

idx_matches            = (d2d_arcsec_closest<posmatch) & (PS_rMag!=-999.) & (PS_iMag!=-999.) & (PS_zMag!=-999.) 

# ATLAS: matches to Pan-STARRS
ATLAS_RA_matches       = ATLAS_RA_closest[idx_matches]
ATLAS_DEC_matches      = ATLAS_DEC_closest[idx_matches]
ATLAS_g_matches        = ATLAS_g_closest[idx_matches]
ATLAS_dg_matches       = ATLAS_dg_closest[idx_matches]
ATLAS_r_matches        = ATLAS_r_closest[idx_matches]
ATLAS_dr_matches       = ATLAS_dr_closest[idx_matches]
ATLAS_i_matches        = ATLAS_i_closest[idx_matches]
ATLAS_di_matches       = ATLAS_di_closest[idx_matches]
ATLAS_z_matches        = ATLAS_z_closest[idx_matches]
ATLAS_dz_matches       = ATLAS_dz_closest[idx_matches]

# Pan-STARRS: matches to ATLA
PS_rMag_matches        = PS_rMag[idx_matches]
PS_rMagErr_matches     = PS_rMagErr[idx_matches]
PS_iMag_matches        = PS_iMag[idx_matches]
PS_iMagErr_matches     = PS_iMagErr[idx_matches]
PS_zMag_matches        = PS_zMag[idx_matches]
PS_zMagErr_matches     = PS_zMagErr[idx_matches]

#########################################################################################################
#                                                                                                       #
#                                           PLOT                                                        #
#                                                                                                       #
#########################################################################################################

# (ATLAS - Pan-STARRS) vs ATLAS

fig = plt.figure(1,figsize=(11,8.5))
ax  = fig.add_subplot(111)
plt.hist2d(ATLAS_r_matches,ATLAS_r_matches-PS_rMag_matches,bins=200,cmin=1)
plt.xlabel("ATLAS r")
plt.ylabel("ATLAS r - Pan-STARRS r")
plt.savefig(plotdir+"ATLAS_vs_PanSTARRS_r.pdf",bbox_inches="tight")
plt.close(fig)

fig = plt.figure(1,figsize=(11,8.5))
ax  = fig.add_subplot(111)
plt.hist2d(ATLAS_i_matches,ATLAS_i_matches-PS_iMag_matches,bins=200,cmin=1)
plt.xlabel("ATLAS i")
plt.ylabel("ATLAS i - Pan-STARRS i")
plt.savefig(plotdir+"ATLAS_vs_PanSTARRS_i.pdf",bbox_inches="tight")
plt.close(fig)

fig = plt.figure(1,figsize=(11,8.5))
ax  = fig.add_subplot(111)
plt.hist2d(ATLAS_z_matches,ATLAS_z_matches-PS_zMag_matches,bins=200,cmin=1)
plt.xlabel("ATLAS z")
plt.ylabel("ATLAS z - Pan-STARRS z")
plt.savefig(plotdir+"ATLAS_vs_PanSTARRS_z.pdf",bbox_inches="tight")
plt.close(fig)


