# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:53:49 2022

@author: Daniel
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def peak_list():
    fish_file_path = 'C:/Users/Daniel/Downloads/Mass_spec_peaks_fish.csv'

    df = pd.read_csv(fish_file_path, header=10)
    df_labels = df.loc[:,['Common name','Reference']]
    df_values = df.iloc[:,4:-2]


    associated_peaks = []
    #Isolate rows from data, remove na values to get list of peaks, sort
    for row in df_values.index:
        #take each row and remove columns with NA values, interpret as list
        peaks = df_values.loc[[row]].dropna(axis=1).values.tolist()
        #sorting values
        sorted_peaks = sorted(peaks)
        #add in species data, sorted_peaks will be list of [peaks, fish name, reference]
        ref = df_labels.loc[row, 'Reference']
        name = df_labels.loc[row, 'Common name']
        sorted_peaks.append(name)
        sorted_peaks.append(ref)
        associated_peaks.append(sorted_peaks)    
    return associated_peaks
    
    