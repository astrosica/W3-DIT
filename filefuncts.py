#!/usr/local/bin/python

#########################################################################################################
#                                                                                                       #
#                                         IMPORT PACKAGES                                               #
#                                                                                                       #
#########################################################################################################


import os
import yaml
import cPickle as pickle


#########################################################################################################
#                                                                                                       #
#                                 IMPORT CONFIGURATION FILE                                             #
#                                                                                                       #
#########################################################################################################


with open("config.yml","r") as fconfig:
    config_data = yaml.safe_load(fconfig)

datadir                  = config_data["datadir"]                  # directory where data is stored
apassdir                 = config_data["apassdir"]                 # directory where APASS catalogues will be stored
dsffdir                  = config_data["dsffdir"]                  # directory where reduced dark-subtracted, flat-fielded science images will be stored
wcsdir                   = config_data["wcsdir"]                   # directory where Astrometry.net WCS solutions will be stored


#########################################################################################################
#                                                                                                       #
#                                       DEFINE FUNCTIONS                                                #
#                                                                                                       #
#########################################################################################################


def direndslash(dir):
    if not dir.endswith("/"):
        dir+="/"
    return dir

def fexists(fpath,f):
    '''
    fname:file to check existence for
    fdir: path to file
    Returns True if file exists, False if it does not
    '''
    if f in os.listdir(fpath):
        return True
    else:
        return False

def pathexists(path):
    '''
    path: directory to check existence for
    Returns True if path exists, False if it does not
    '''
    if os.path.exists(path):
        return True
    else:
        return False

def bothcaldirs(caldir,flatdir):
        if not (
            pathexists(caldir) and pathexists(flatdir)
            ):
            if not pathexists(caldir) and not pathexists(flatdir):# if both missing
                return False
            elif not pathexists(caldir) and pathexists(flatdir): # if only calibration folder missing
                return "Calibration"
            elif pathexists(caldir) and not pathexists(flatdir): # if only flat folder missing
                return "AutoFlat"
        else:
            return True

def fcheckext(file,folder,ext,wrongextfiles):
    if file.endswith(ext):
        return True, wrongextfiles
    else:
        wrongextfiles = appenditemtodict(file,folder,wrongextfiles)
        return False, wrongextfiles

def fcheckflag(file,folder,badfiles,logstring,store,log=False):
    filesplit = re.split("-|\.",file)
    if FLAG and flag in filesplit:
        if store:
            badfiles = appenditemtodict(file,folder,badfiles)
        return True, badfiles
    else:
        return False, badfiles

def fcheckobj(file,folder,badobjfiles,logstring,store,log=False):
    filesplit = re.split("-|\.",file)
    temp = 0
    for object in objects:
        if object in filesplit:
            temp+=1
    if temp==0:
        if file!=calfolder[:-1] and file!=flatfolder[:-1] and store==True:
            badobjfiles = appenditemtodict(file,folder,badobjfiles)
        return False,badobjfiles
    elif temp==1:
        return True,badobjfiles
    else:
        if log:
            logstring.append("    > WARNING: more than one object in filename")

def findfolderdir(folder,VERBOSE=False):
    apassfolderdir  = apassdir+folder+"/" # folder to contain APASS catalogues
    dsfffolderdir   = dsffdir+folder+"/"  # folder containing dark subtracted flat fielded science images
    wcsfolderdir    = wcsdir+folder+"/"   # folder containing original wcs-solved science images
    datafolderdir   = datadir+folder+"/"  # folder containing original raw science images
    if pathexists(dsfffolderdir):
        folderdir = dsfffolderdir
        if VERBOSE:
            print "processing dark-subtracted flat-fielded science images"
    elif pathexists(wcsfolderdir):
        folderdir = wcsfolderdir
        if VERBOSE:
            print "processing original WCS-solved science images"
    elif pathexists(datafolderdir):
        folderdir = datafolderdir
        if VERBOSE:
            print "processing original raw science images"
    
    return folderdir

def nelementsindict(dict):
    nelements = []
    keys = dict.keys()
    for key in keys:
        elements = dict[key]
        n = len(elements)
        nelements.append(n)
    return nelements

def fwritedict_pickle(dict,filedir):
    '''
    Writes the contents of a dictioinary to a pickled file.
    '''

    with open(filedir,"w+") as file:
            file.write(pickle.dumps(dict))

def freaddict_pickle(filedir):
    '''
    Reads the contents of a dictioinary from a pickled file.
    '''

    with open(filedir,"rb") as infile:
        dict = pickle.load(infile)

    return dict

def fwritedict_text(dict,filedir):
    '''
    Writes the contents of a dictioinary to a regular text file.
    '''

    with open(filedir,"w+") as f:
        for key in dict.keys():
            f.writelines(key+"\n")
            filelist = dict[key]
            for file in filelist:
                f.writelines(file+"\n")


