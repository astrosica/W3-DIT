#!/usr/local/bin/python

import os
import re
import yaml
import numpy as np
from astropy.io import fits
from astropy.io.fits import getdata

import imp
# filefuncts
import filefuncts
imp.reload(filefuncts)
from filefuncts import *
# collectzipdata
import collectzipdata
imp.reload(collectzipdata)
from collectzipdata import *

import multiprocessing
from functools import partial

########## IMPORT CONFIG FILE ##########

with open("config.yml","r") as fconfig:
    config_data = yaml.safe_load(fconfig)

ncores_cal               = int(config_data["ncores_cal"])
datadir                  = config_data["datadir"]
dsffdir                  = config_data["dsffdir"]
mastercaldir             = config_data["mastercaldir"]
MMcaldir                 = mastercaldir+"MMcalimages/"
sciencemaskdir           = config_data["sciencemaskdir"]
FLAG                     = config_data["FLAG"]
flag                     = config_data["flag"]
objects                  = config_data["objects"].split(" ")
calfolder                = config_data["calfolder"]
flatfolder               = config_data["flatfolder"]
sciencedata_dsff_cutoff  = config_data["sciencedata_dsff_cutoff"]
wcsdir                   = config_data["wcsdir"]

if sciencedata_dsff_cutoff!=None:
    sciencedata_dsff_cutoff = float(sciencedata_dsff_cutoff)

pool_cal = multiprocessing.Pool(ncores_cal)

############## FUNCTIONS ##############

def createmasterbias(caldir,folder,badbiasfiles,logstring,store,log=False):
    allbias = []
    for file in os.listdir(caldir):
        if file.startswith("Bias"):
            checkflag,badbiasfiles=fcheckflag(file,folder,badbiasfiles,logstring,store,log=False)
            if checkflag:
                if log:
                    logstring.append("    > skipping bad file: "+file)
            else:
                bias=getdata(caldir+file)
                allbias.append(bias)
    masterbias = np.median(np.dstack(allbias),axis=2)
    return masterbias,badbiasfiles

def createMMbias(all_masterbias):
    MMbias = np.median(np.dstack(all_masterbias),axis=2)
    return MMbias

def createmasterdarks(caldir,masterbias,folder,baddarkfiles,logstring,store,log=False):
    darkdict = {}
    masterdarks = {}
    for file in os.listdir(caldir):
        if file.startswith("Dark"):
            checkflag,baddarkfiles=fcheckflag(file,folder,baddarkfiles,logstring,store,log=False)
            if checkflag:
                if log:
                    logstring.append("    > skipping bad file: "+file)
            else:
                dark,darkheader=getdata(caldir+file,header=True)
                dark_bs = dark - masterbias
                exptime=darkheader["EXPTIME"]
                darkdict = appenditemtodict(dark_bs,exptime,darkdict)
    for darkexptime in darkdict.keys():
        darks = darkdict[darkexptime]
        masterdark = np.median(np.dstack(darks),axis=2)
        masterdarks[darkexptime] = masterdark
    return masterdarks,baddarkfiles

def appendmasterdarks(masterdarks,all_masterdarks):
    for darkkey in masterdarks.keys():
        masterdark = masterdarks[darkkey]
        all_masterdarks = appenditemtodict(masterdark,darkkey,all_masterdarks)
    return all_masterdarks

def createMMdarks(all_masterdarks):
    MMdarks = {}
    for darkexptime in all_masterdarks.keys():
        masterdarks=all_masterdarks[darkexptime]
        MMdark = np.median(np.dstack(masterdarks),axis=2)
        MMdarks[darkexptime]=MMdark
    return MMdarks

def missingcalibration(folder,folderdir,badsciencefiles,badobjfiles,missing_calibration,logstring,store,log=False):
    for file in os.listdir(folderdir):
        checkobj,badobjfiles=fcheckobj(file,folder,badobjfiles,logstring,store,log=False)
        if checkobj==True:
            checkflag,badsciencefiles = fcheckflag(file,folder,badsciencefiles,logstring,store,log=False)
            missing_calibration = appenditemtodict(file,folder,missing_calibration)

def findclosestexptime(exptime,exptimes):
    absdiff = np.abs(np.array(exptimes)-exptime)
    wheremin = np.where(absdiff==absdiff.min())[0][0]
    closestexptime = exptimes[wheremin]
    return closestexptime

def createmasterflats(flatdir,masterbias,masterdarks,folder,badflatfiles,logstring,store,log=False):
    darkexptimes = masterdarks.keys()
    darkexptimes.sort()
    darkexptimes.insert(0,0.0)
    masterflats = {}
    flatdict = {}
    for file in os.listdir(flatdir):
        if file.startswith("AutoFlat"):
            checkflag,badflatfiles=fcheckflag(file,folder,badflatfiles,logstring,store,log=False)
            if checkflag:
                if log:
                    logstring.append("    > skipping bad file: "+file)
            else:
                flatdata,flatheader=getdata(flatdir+file,header=True)
                flatexptime=flatheader["EXPTIME"]
                passband=flatheader["FILTER"]
                darkexptime = findclosestexptime(flatexptime,darkexptimes)
                if darkexptime<0.00001:        # requires only bias subtraction
                    darksubflat = flatdata-masterbias
                else:                          # requires both bias and dark subtraction
                    masterdark = masterdarks[darkexptime]
                    darksubflat = flatdata-masterbias-masterdark
                darksubflat_norm = darksubflat/np.median(darksubflat)
                flatdict = appenditemtodict(darksubflat_norm,passband,flatdict)
    for passband in flatdict.keys():
        flats = flatdict[passband]
        masterflat = np.median(np.dstack(flats),axis=2)
        masterflat_norm = masterflat/np.median(masterflat)
        masterflats[passband] = masterflat_norm
    return masterflats,badflatfiles

def appendmasterflats(masterflats,all_masterflats):
    for flatkey in masterflats.keys():
        masterflat = masterflats[flatkey]
        all_masterflats = appenditemtodict(masterflat,flatkey,all_masterflats)
    return all_masterflats

def createMMflats(all_masterflats):
    MMflats = {}
    for passband in all_masterflats.keys():
        masterflats=all_masterflats[passband]
        MMflat = np.median(np.dstack(masterflats),axis=2)
        MMflat_norm = MMflat/np.median(MMflat)
        MMflats[passband]=MMflat_norm
    return MMflats

def missingflat(folder,folderdir,badsciencefiles,badobjfiles,missing_flat,logstring,store,log=False):
    for file in os.listdir(folderdir):
        checkobj,badobjfiles=fcheckobj(file,folder,badobjfiles,logstring,store,log=False)
        if checkobj==True:
            checkflag,badsciencefiles = fcheckflag(file,folder,badsciencefiles,logstring,store,log=False)
            missing_flat = appenditemtodict(file,folder,missing_flat)

def missingcalflat(folder,folderdir,badsciencefiles,badobjfiles,missing_calflat,logstring,store,log=False):
    for file in os.listdir(folderdir):
        checkobj,badobjfiles=fcheckobj(file,folder,badobjfiles,logstring,store,log=False)
        if checkobj==True:
            checkflag,badsciencefiles = fcheckflag(file,folder,badsciencefiles,logstring,store,log=False)
            missing_calflat = appenditemtodict(file,folder,missing_calflat)

def savemasterbias(mastercalfolderdir,masterbias):
    fits.writeto(mastercalfolderdir+"masterbias.fts",masterbias,clobber=True)

def savemasterdarks(mastercalfolderdir,masterdarks):
    for darkexptime in masterdarks.keys():
        masterdark = masterdarks[darkexptime]
        fits.writeto(mastercalfolderdir+"masterdark-"+str(int(darkexptime))+".fts",masterdark,clobber=True)

def savemasterflats(mastercalfolderdir,masterflats):
    for passband in masterflats.keys():
        masterflat = masterflats[passband]
        fits.writeto(mastercalfolderdir+"masterflat-"+passband+".fts",masterflat,clobber=True)

def bdsff(data,exptime,passband,masterbias,masterdarks,masterflats):
    # bias/dark subtraction
    darkexptimes = masterdarks.keys()
    darkexptimes.sort()
    darkexptimes.insert(0,0.0)
    darkexptime = findclosestexptime(exptime,darkexptimes)
    if darkexptime<0.00001:
        data_bds = data-masterbias
    else:
        data_bs = data-masterbias
        masterdark = masterdarks[darkexptime]
        data_bds = data_bs-masterdark
    # flat fielding
    masterflat = masterflats[passband]
    data_bdsff = data_bds/masterflat
    return data_bdsff

def MMbdsff(data,mastercalfolderdir,exptime,passband,mastermastercaldir,useMMbias=False,useMMdark=False,useMMflat=False):
    # bias
    if useMMbias:
        masterbias=getdata(mastermastercaldir+"masterbias.fts")
        data_bs = data-masterbias
    else:
        masterbias=getdata(mastercalfolderdir+"masterbias.fts")
    # dark
    if useMMdark:
        masterdark=getdata(mastermastercaldir+"masterdark-"+str(int(exptime))+".fts")
        data_bds = data_bs-masterdark
    else:
        darkexptimes = []
        for file in os.listdir(mastercalfolderdir):
            if file.startswith("masterdark"):
                filesplit = re.split("-|\.",file)
                darkexptime = float(filesplit[-2])
                darkexptimes.append(darkexptime)
        darkexptimes.sort()
        darkexptimes.insert(0,0.0)
        darkexptime = findclosestexptime(exptime,darkexptimes)
        if darkexptime<0.00001:
            data_bds = data-masterbias
        else:
            data_bs = data-masterbias
            masterdark = getdata(mastercalfolderdir+"masterdark-"+str(int(darkexptime))+".fts")
            data_bds = data_bs-masterdark
    # flat fielding
    if useMMflat:
        masterflat=getdata(mastermastercaldir+"masterflat-"+passband+".fts")
    else:
        masterflat=getdata(mastercalfolderdir+"masterflat-"+passband+".fts")
    data_bdsff = data_bds/masterflat
    return data_bdsff

def bdsff_missingdata(
    missing_dict,VERBOSE,GENERATE,logstring,badsciencefiles,badobjfiles,wrongextfiles,log=False,useMMbias=False,useMMdark=False,useMMflat=False
    ):
    for folder in missing_dict.keys():
        if VERBOSE==True:
            print "processing folder "+folder+"..."
        files = missing_dict[folder]
        if len(files)>=ncores_cal:
            bdsff_missingdata_f_mp = partial(bdsff_missingdata_f,folder,VERBOSE,GENERATE,log=log,useMMbias=useMMbias,useMMdark=useMMdark,useMMflat=useMMflat)
            badsciencefiles_temp_zip,badobjfiles_temp_zip,wrongextfiles_temp_zip,logstring_temp_zip=zip(*pool_cal.map(bdsff_missingdata_f_mp,files))
            # collect zipped lists and dictionaries
            badsciencefiles_temp=collectziplists(badsciencefiles_temp_zip)
            badobjfiles_temp=collectziplists(badobjfiles_temp_zip)
            wrongextfiles_temp=collectzipdicts(wrongextfiles_temp_zip)
            logstring_temp=collectziplists(logstring_temp_zip)
            # append
            for tempfile in badsciencefiles_temp:
                badsciencefiles.append(tempfile)
            for tempfile in badobjfiles_temp:
                badobjfiles.append(tempfile)
            for tempkey in wrongextfiles_temp:
                if tempkey not in wrongextfiles.keys():
                    wrongextfiles.setdefault(tempkey,[])
                tempfiles = wrongextfiles_temp[tempkey]
                for tempfile in tempfiles:
                    wrongextfiles[tempkey].append(tempfile)
        else:
            for file in files:
                badsciencefiles_temp,badobjfiles_temp,wrongextfiles_temp,logstring_temp=bdsff_missingdata_f(
                    folder,VERBOSE,GENERATE,file,log,useMMbias,useMMdark,useMMflat
                    )
            # append
            for tempfile in badsciencefiles_temp:
                badsciencefiles.append(tempfile)
            for tempfile in badobjfiles_temp:
                badobjfiles.append(tempfile)
            for tempkey in wrongextfiles_temp:
                if tempkey not in wrongextfiles.keys():
                    wrongextfiles.setdefault(tempkey,[])
                tempfiles = wrongextfiles_temp[tempkey]
                for tempfile in tempfiles:
                    wrongextfiles[tempkey].append(tempfile)
    return badsciencefiles,badobjfiles,wrongextfiles,logstring

def bdsff_missingdata_f(folder,VERBOSE,GENERATE,file,log,useMMbias,useMMdark,useMMflat):
    logstring_temp = []
    badsciencefiles_temp = []
    badobjfiles_temp = []
    wrongextfiles_temp = {}
    wcsfolderdir=wcsdir+folder+"/"
    if pathexists(wcsfolderdir):
        folderdir=wcsdir+folder+"/"
    else:
        folderdir = datadir+folder+"/"                   # folder containing original raw data
    dsfffolderdir = dsffdir+folder+"/"                   # folder containing dark subtracted flat fielded science images
    mastercalfolderdir = mastercaldir+folder+"/"         # folder containing master calibration images
    sciencemaskfolderdir = sciencemaskdir+folder+"/"     # folder containing science masks
    #
    checkobj,badobjfiles_temp=fcheckobj(file,folder,badobjfiles_temp,logstring_temp,store=True,log=log)
    checkext,wrongextfiles_temp=fcheckext(file,folder,".fts",wrongextfiles_temp)
    if checkobj and checkext:
        checkflag,badsciencefiles_temp=fcheckflag(file,folder,badsciencefiles_temp,logstring_temp,store=True,log=log)
        if checkflag:
            if log:
                logstring_temp.append("    > skipping bad file: "+file)
            if VERBOSE:
                print "    > skipping bad file: "+file
        else:
            checkobj,badobjfiles_temp=fcheckobj(file,folder,badobjfiles_temp,logstring_temp,store=True,log=log)
            if checkobj:
                if VERBOSE:
                    print "    > processing file: "+file
                if log:
                    logstring_temp.append("    > processing file: "+file)
    #
                sciencedata, scienceheader=getdata(folderdir+file,header=True)
                exptime=scienceheader["EXPTIME"]
                passband=scienceheader["FILTER"]
                sciencedata_bdsff = MMbdsff(sciencedata,mastercalfolderdir,exptime,passband,MMcaldir,useMMbias=useMMbias,useMMdark=useMMdark,useMMflat=useMMflat)
                if sciencedata_dsff_cutoff!=None and sciencedata_bdsff.min()<sciencedata_dsff_cutoff:
                    if log:
                        logstring_temp.append(file+ "was masked")
                    scienceheader["MASK"]="True"
                    scienceheader["MASKVAL"]=sciencedata_dsff_cutoff
                    sciencedata_bdsff,sciencemask_bdsff = maskdata(sciencedata_bdsff,sciencedata_dsff_cutoff)
                    if GENERATE:
                        savemask(sciencemaskfolderdir,file,sciencemask_bdsff)
                savebdsff(dsfffolderdir,file,sciencedata_bdsff,scienceheader)
    return badsciencefiles_temp,badobjfiles_temp,wrongextfiles_temp,logstring_temp

def savebdsff(dsfffolderdir,file,sciencedata_bdsff,scienceheader):
    newfile = file.replace(".fts","-dsff.fts")
    fits.writeto(dsfffolderdir+newfile,sciencedata_bdsff,scienceheader,clobber=True)

def maskdata(data,cutoff):
    above = np.where(data>=cutoff)
    below = np.where(data<cutoff)
    mask = np.copy(data)
    mask[above] = 1.0
    mask[below] = 0.0
    maskeddata = data*mask
    return [maskeddata, mask]

def savemask(maskdir,file,mask):
    newfile = file.replace(".fts","-dsff-mask.fts")
    fits.writeto(maskdir+newfile,mask,clobber=True)
