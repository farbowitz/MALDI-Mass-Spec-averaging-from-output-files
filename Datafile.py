# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 20:51:54 2022

@author: Daniel
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import argrelmax

#ATTEMPTING LOGGING
import logging
import sys

logger = logging.getLogger(__name__)
stream_handler = logging.StreamHandler(sys.stdout)
logger.addHandler(stream_handler)

known_peaks_path="C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Mass_spec_peaks_fish.csv"
project_path = "C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Usable Spectra July 2022/"

def get_all_files(path):
    file_list = []
    for (dirpath, dirnames, filenames) in os.walk(path):
        file_list.extend(filenames)
    return file_list

##find all text files within current directory
def get_txt_files(path):
    file_list = []
    for element in os.listdir(path):
        if element.endswith('.txt'):
            file_list.append(element)
    return file_list

def list_similarity_operation(list1, list2, allowable_distance=0.5):
    match_list = []
    non_matches1 = []
    non_matches2 = []
    if len(list1) ==0 or len(list2) == 0:
        print('No peaks to form comparison from. :(')
        return None, None, None
    for x_value in list1:
        if find_nearest(list2, x_value, None) < allowable_distance:
            index1 = np.where(list1 == x_value)[0]
            index2 = find_nearest_index(list2, x_value)
            match_list.append(np.around(float(list1[index1]),2))#, float(list2[index2])))
            non_matches2 = np.delete(list2, index2)
        else:
            non_matches1.append(x_value)
    return match_list, non_matches1, non_matches2

#find_nearest taken from TA program, useful for following function
def find_nearest(array, number, direction): 
    if direction is None:
        idx = (np.abs(array - number)).min()
    elif direction == 'backward':
        _delta = number - array
        _delta_positive = _delta[_delta > 0]
        if not _delta_positive.empty:
            idx = _delta_positive.min()
    elif direction == 'forward':
        _delta = array - number
        _delta_positive = _delta[_delta >= 0]
        if not _delta_positive.empty:
            idx = _delta_positive.min()
    return idx

def find_nearest_index(array, number):
    delta_array = np.abs(array-number)
    index = np.where(delta_array == delta_array.min())[0]
    return index


## individual txt files from each run are a datafile
#Which direction should inheritance run
class Datafile:
    def __init__(self, txt_filepath):
        #import data into pandas based on filename referenced
        self.filepath = txt_filepath

        df = pd.read_csv(txt_filepath, sep="\t", header=0)
        df.columns = ['m/Z', 'intensity']
        self.DataFrame = self.clip_data(df)
        #return none if no data
        if not any(self.DataFrame):
            logger.debug('No data available in dataframe.')
        self.props_from_filepath()

        self.DataFrame



    #PROVIDES RUN_NAME, CELL, SAMPLE_ID    
    def props_from_filepath(self):
        path_split = self.filepath.split(sep='/')
        run_folder = path_split[-2]
        filename = path_split[-1]
        self.run_name = os.path.commonprefix([run_folder, filename])
        #need to remove initial date
        self.run_name = '_'.join(self.run_name.split(sep='_')[1:])
        #if not self.run_name:
            #logger.warning('No common name between file and folder. Check {}'.format(run_folder))
        self.sample_id = filename.removesuffix('.txt')
        
        #for other filetypes
        #self.date = file_comps[0]
        #self.cell = file_comps[-1]
        #check cell matches format with regex?
        #self.sample_id = file_comps[-2:]
        #logger.info('Creating Dataframe for cell {} (sample {}) in {}'.format(self.cell, self.sample_id, self.run_name))

    def locate_run_map(self, folder='maps'):
        map_folder = project_path+folder+'/'
        if os.path.exists(map_folder):
            #if filename matches run name
            maps = get_all_files(map_folder)
            for map in maps:
                #better method needed, not index based (use split and check if lists are equivalent?), removed extra underscore from end of name
                if self.run_name in map:
                  map_path = map_folder+map #filename, once found
            try:
                return map_path
            except:
                logger.warning('No map found containing {} in {}.'.format(self.run_name, map_folder))
            
        else:
            logger.warning('No maps folder found. Creating one at {}'.format(map_folder))
            #create folder
            #add default.csv file to it
            map_path=None
        return map_path

    def import_run_map(self, folder='maps'):
        map_path = self.locate_run_map()
        df = pd.read_csv(map_path, names=['cell id', 'sample id'], header=None)
        #df is assumed to have the first column be cell number, second column be sample name (CAL for calibrant)
        df = df.dropna()
        df = df.drop(df[df['sample id'] == 'CAL'].index)
        #remove nan values, ignore CAL
        return df

    def get_sample_id(self):
        #load pandas dataframe associated with this run
        df = self.import_run_map()
        #use .item() in pandas to return single object as original datatype rather than pandas object
        sample_id = df.loc[df['cell id']==self.cell, 'sample id'].item()
        return str(sample_id)

    #DATA RELATED FUNCTIONS
    #select region of interest, for now default is (850, 3000)
    def clip_data(self, df, min=850, max=3000):
        clipped_df = df[df['m/Z'].between(min, max)]
        return clipped_df

    def normalize(self):
        pass

    ##PEAK FINDING
    
    def identify_peaks(self):
        return None

    def peak_coords(self):
        df = self.DataFrame
        #need numpy array to use argrelmax?
        y = np.asarray( df['intensity'])
        x = np.asarray( df['m/Z'])

        #prominence value based on max for time used, or over whole dataset?
        int_max = y.max() 
        #normalize
        y = y/int_max
        #peaks, other = scipy.signal.find_peaks(y, prominence=0.05*df_max)
        #alternatively, use argrelextrema?
        n=200
        indices = argrelmax(y, order=n)
        x_peaks = x[indices]
        y_peaks = y[indices]
        return x_peaks, y_peaks

    def plot_peaks(self):
        x_peaks, y_peaks = self.peak_coords()
        #may want to get use out of x_peaks info
        plt.plot(x_peaks, y_peaks, 'xb')
        #add in labels?
        
        for peak in x_peaks:
            plt.annotate(round(peak), xy=(peak, y_peaks[np.where(x_peaks == peak)]))
        #plt.vlines(x_peaks, ymin=0, ymax=y_peaks)
        self.plot()

    def list_peaks(self):
        x_peaks, y_peaks = self.peak_coords()
        print(x_peaks)
        return None

    def get_peak_list(self):
        x_peaks, y_peaks = self.peak_coords()
        return x_peaks
    
    def clean_noise(self):
        return None

    #internal function to plot data, without_show exists for comparative purposes
    def plot_without_show(self):
        plt.plot(self.DataFrame['m/Z'], self.DataFrame['intensity']/self.DataFrame['intensity'].max())
        plt.xlabel('m/Z')
        plt.ylabel('Intensity (a.u.)')
    
    def plot(self):
        self.plot_without_show()
        plt.show()



    def load_known_peaks(self, known_peaks_path="C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Mass_spec_peaks_fish.csv"):
        df = pd.read_csv(known_peaks_path, header=10)
        df_vals = df[df.columns[4:-2]]
        df_meta = df[df.columns[[0,1,2,3,-2,-1]]]

        #indices should match between list and df_meta to correctly assign fish names
        list_of_known_peaks = [np.asarray(df_vals.loc[i].dropna()) for i in range(len(df_vals))]
        return list_of_known_peaks, df_meta
    
    def list_matches_to_fish_data(self,file_peaks):
        list_of_matches =[]
        list_of_known_peaks, df_meta = self.load_known_peaks()
        for known_peaks in list_of_known_peaks:
            matches, remainder1, remainder2 = list_similarity_operation(file_peaks, known_peaks)
            list_of_matches.append(matches)
        match_lengths = [len(match_list) for match_list in list_of_matches]
        return list_of_matches, match_lengths, df_meta

    def print_species_matches_from_datafile(self,path):
        obj = Datafile(path)
        file_peaks = obj.get_peak_list()
        print('Sample: {}'.format(obj.sample_id))
        list_of_matches, match_lengths, df_meta = self.list_matches_to_fish_data(file_peaks)
        #list the length of each matches list to align indices

        #assuming this puts indices in reverse order of highest
        indices = np.argpartition(match_lengths, -3)[-3:]
        print('Highest number of peak matches:')
        for index in indices:
            top_val = 3-np.where(indices==index)[0]
            top_spec = df_meta.loc[index]['Species']
            top_len = match_lengths[index]
            print('#{}: {}, {} matches @ {}'.format(top_val, top_spec, top_len, list_of_matches[index]))
        logger.info('Datafile analysis complete. \n')
        

        
    
'''
#test code for datafile class
test_path = 'C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Liz_MALDI_runs/20220711_LMQ_Test Text Files/20220711_LMQ_Test_300622_A3.txt'
test_path=test_path.replace('A3', 'C2')
obj = Datafile(test_path)
obj.print_species_matches_from_datafile(test_path)
'''