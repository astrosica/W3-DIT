#!/usr/local/bin/python

import numpy as np

def collectziplists(ziplist):
    list = []
    for listitem in ziplist:
        if listitem!=[]:
            for i in range(len(listitem)):
                item = listitem[i]
                list.append(item)
    return list

def collectzipdicts(zipdict):
    dictionary = {}
    for dict in zipdict:
        keys = dict.keys()
        for key in keys:
            if key not in dictionary.keys():
                dictionary.setdefault(key,[])
            for i in range(len(dict[key])):
                data=dict[key][i]
                dictionary[key].append(data)
    return dictionary

def appenditemtodict(item,key,dict):
    if key not in dict.keys():
        dict.setdefault(key,[])
    dict[key].append(item)
    
    return dict

def append_dict2dict(dict1,dict2):
    '''
    Appends contents of dict2 to dict1.
    '''
    
    for dict2_key in dict2.keys():
        if dict2_key not in dict1.keys():
            dict1.setdefault(dict2_key,[])
        dict1[dict2_key].append(dict2[dict2_key][0]) # modifier
    
    return dict1

def append_dict2dict_temp(dict1,dict2):
    '''
    Appends contents of dict2 to dict1.
    '''
    
    for dict2_key in dict2.keys():
        if dict2_key not in dict1.keys():
            dict1.setdefault(dict2_key,[])
        dict1[dict2_key].append(dict2[dict2_key])
    
    return dict1

def flattendict(dict,ARRAYS=False):
    dict_flat = {}
    for key in dict.keys():
        data = dict[key]
        if ARRAYS==True:
            data_flat = np.array([item for sublist in data for item in sublist])
        else:
            data_flat = np.array(data)
        dict_flat.setdefault(key,data_flat)
    
    return dict_flat


