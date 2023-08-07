#!/usr/local/bin/python

'''
callAPASS - contains a function to create APASS catalogues
Requires the following modules: urllib, docopt, numpy, os
Contains the following functions: callAPASS

Originally adapted from Natalie Price-Jones
'''


########################## IMPORT PACKAGES ###########################


import os
import numpy as np
import urllib as url

import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


########################## FUNCTIONS ###########################

def callapass(file,ra,dec,fov,apassdir,VERBOSE,GENERATE,totalfovbool,useobjpos=True):
    '''
    Finds an appropriate APASS catalogue for an image and saves it in apassdir
    ra:         right ascension in degrees
    dec:        declination in degrees
    fov:        field of view in degrees (max 15)
    apassdir:   directory in which apass files are located
    
    Saves an array (apassdata), where each column contains a column of APASS data
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
    
    # Create URL by which to access catalogue
    address="http://www.aavso.org/cgi-bin/apass_download.pl?ra={0}&dec={1}&radius={2}&outtype=0".format(ra,dec,fov/2.0)
    catalog = url.urlopen(address)
    
    # Do a lot of terrible data parsing
    apassdata = catalog.read()
    catalog.close()
    apassdata = apassdata.replace("<Td><font size=-1>","")
    apassdata = apassdata.replace("</td>","")
    apassdata = apassdata[577:len(apassdata)-256]
    apassdata = apassdata.split("\n<tr>")
    apassdata = np.array(apassdata)[np.where(np.array(apassdata) != "")]
    apassdata = np.array([i.replace(" ","").replace("NA","-1").split("\n\t") for i in apassdata])
    apassdata = np.delete(apassdata,5,axis=1)
    apassdata = apassdata.astype("float")
    
    # Save catalogue to file
    if totalfovbool==True:
        apassfname = apassdir+"ra{0}dec{1}fov{2}.apass".format(ra,dec,fov)
    else:
        apassfname = apassdir+file.replace(".fts",".apass")
    if VERBOSE:
        print "creating APASS catalogue: "+apassfname
    np.savetxt(apassfname,apassdata)

    if GENERATE:

        if VERBOSE:
            print "creating annotation files..."

        # COVERAGE

        color_kvis = "RED"
        color_ds9 = color_kvis.lower()
        size = fov/2.
        fname_kvis = apassfname.replace(".apass","-coverage.ann")
        fname_ds9 = apassfname.replace(".apass","-coverage.reg")
        f_kvis = open(fname_kvis,"w+")
        f_ds9 = open(fname_ds9,"w+")
        f_kvis.write("COLOR %s\n\n" % color_kvis)
        f_kvis.write("CIRCLE {0} {1} {2}\n\n".format(ra,dec,size))
        f_ds9.write("fk5;circle {0} {1} {2} # color={3}\n\n".format(ra,dec,size,color_ds9))
        f_kvis.close()
        f_ds9.close()

        if VERBOSE:
            print "wrote kvis annotation file "+fname_kvis
            print "wrote ds9 annotation file "+fname_ds9

        # CATALOGUE

        ra_all = apassdata[:,0]
        ra_err_all = apassdata[:,1]/3600. # arcsec --> deg
        dec_all = apassdata[:,2]
        dec_err_all = apassdata[:,3]/3600. # arcsec -->deg

        color_kvis = "RED"
        color_ds9 = color_kvis.lower()
        fname_kvis = apassfname.replace(".apass",".ann")
        fname_ds9 = apassfname.replace(".apass",".reg")
        f_kvis = open(fname_kvis,"w+")
        f_ds9 = open(fname_ds9,"w+")
        f_kvis.write("COLOR %s\n\n" % color_kvis)
        for i in range(len(ra_all)):
            ra_i = ra_all[i]
            ra_err_i = ra_err_all[i]
            dec_i = dec_all[i]
            dec_err_i = dec_err_all[i]
            pos_err = np.sqrt((ra_err_i)**2. + (dec_err_i)**2.)
            pos_err = 0.002
            f_kvis.write("CIRCLE {0} {1} {2}\n\n".format(ra_i,dec_i,pos_err))
            f_ds9.write("fk5;circle {0} {1} {2} # color={3}\n\n".format(ra_i,dec_i,pos_err,color_ds9))
        f_kvis.close()
        f_ds9.close()

        if VERBOSE:
            print "wrote kvis annotation file "+fname_kvis
            print "wrote ds9 annotation file "+fname_ds9
    
    return apassdata
