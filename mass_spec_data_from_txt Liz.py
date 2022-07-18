# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 19:15:28 2022

@author: Daniel
"""
# LIBRARY DEPENDENCIES

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

# FUNDAMENTAL VARIABLE INPUTS
known_peaks_path="C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Mass_spec_peaks_fish.csv"
project_path = "C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Liz_MALDI_runs/"
## OPTIONAL- SAMPLE NAMING METHOD: assign a folder to keep CSV files (titled with the run name) which map column 1 wells (e.g. A1, B3, etc.) to sample names (e.g. sample 2, coppergate_fish). Use CAL for calibrant wells. 
maps_folder_name = 'maps'
## default - adds folder for sample data under projet directory

###TEST DATA
test_folder = project_path+'20220429_Coppergate_1 Text Files/'
test_filename = '20220429_Coppergate_1_A8.txt'
test_path = test_folder+test_filename

#USEFUL FUNCTIONS

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

# CLASSES

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
        if not self.run_name:
            logger.warning('No common name between file and folder. Check {}'.format(run_folder))
        filename = filename.removesuffix('.txt')
        file_comps =filename.split(sep='_')
        self.date = file_comps[0]
        self.cell = file_comps[-1]
        #check cell matches format with regex?
        self.sample_id = self.get_sample_id()
        logger.info('Creating Dataframe for cell {} (sample {}) in {}'.format(self.cell, self.sample_id, self.run_name))

    def locate_run_map(self, folder=maps_folder_name):
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
                logger.warning('No map found containing {} in {}.'.format(self.run_name, ))
            
        else:
            logger.warning('No maps folder found. Creating one at {}'.format(map_folder))
            #create folder
            #add default.csv file to it
            map_path=None
        return map_path

    def import_run_map(self, folder=maps_folder_name):
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
        plt.plot(self.DataFrame['m/Z'], self.DataFrame['intensity'])
        plt.xlabel('m/Z')
        plt.ylabel('Intensity (a.u.)')
    
    def plot(self):
        self.plot_without_show()
        plt.show()


#test code for datafile class
#test_path=test_path.replace('A2', 'C2')
obj = Datafile(test_path)
print(obj.cell)
obj.plot_peaks()

#New class for peak lists?

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


#put into run, but try to think of appropriate metric
#maybe do pairwise similar elements where you pop them out of their respective lists and see what remains
def list_similarity_metric(list_a, list_b):
    unique_elements = 0
    dev = 0
    for x_value in list_a:
        diff = find_nearest(list_b, x_value, direction=None)
        sq_diff = diff**2
        #assign unique element for 
        if sq_diff > 4:
            unique_elements += 1
        else:
            dev += sq_diff
    pct_unique = 100*unique_elements/len(list_a)
    mtr = 1-(sq_diff/((len(list_a)-unique_elements)*4))
    print('The first list has {} percent unique elements, similar elements have a metric of {}'.format(pct_unique,mtr))

def list_similarity_operation(list1, list2):
    match_list = []
    new_list1 = []
    if len(list1) ==0 or len(list2) == 0:
        print('No peaks to form comparison from. :(')
        return None, None, None
    for x_value in list1:
        if find_nearest(list2, x_value, None) < 1:
            index1 = np.where(list1 == x_value)[0]
            index2 = find_nearest_index(list2, x_value)
            match_list.append((float(list1[index1])))#, float(list2[index2])))
            list2 = np.delete(list2, index2)
        else:
            new_list1.append(x_value)
    return match_list, new_list1, list2



'''
test_path_2 = test_path.replace('C2', 'A2')

obj2 = Datafile(test_path_2)
obj2.list_peaks()
list_a = obj.get_peak_list()
list_b = obj2.get_peak_list()
list_similarity_metric(list_a, list_b)
match, list1, list2 = list_similarity_operation(list_a, list_b)
print(list1)
print(list2)
print(match)
print(len(list_a))
print(len(list_b))
print(len(list1))
print(len(list2))
print(len(match))
print([np.abs(el1-el2) for (el1, el2) in match])
'''



## Creates a Sample from one or multiple datafiles
class Sample:
    def __init__(self, list_of_datafiles, output_name, folder_path, output_path):
        self.size = len(list_of_datafiles)
        self.list = list_of_datafiles
        self.name = output_name
        self.folder_path = folder_path
        self.output_path = output_path
        #make empty file of column length, *SHOULD ASSIGN USING ARRAY LENGTH, UNLESS ALL MALDI RUNS ARE SAME SIZE* 
        running_total = [0] * 107875
        run_file_list = get_txt_files(self.folder_path)
        for file in run_file_list:
            #check if dataset file has same length and index
            filename = Title(file)
            if (str(filename.row)+str(filename.column)) in self.list:
               running_total += np.asarray(Datafile(self.folder_path+file).DataFrame['intensity'])
        df_average = running_total/self.size
        df_index = Datafile(self.folder_path+file).DataFrame['m/Z']
        df = pd.DataFrame(df_average, index=df_index)
        self.DataFrame = df
        
    def output_to_file(self):
        self.DataFrame.to_csv(r''+self.output_path+self.name+'.txt', sep='\t', index=True, header=False)
        #check if liz wants headers removed
        
    def identify_peaks(self):
        return None
    
    def clean_noise(self):
        return None        
        
    def plot(self):
        plt.plot(self.DataFrame['m/Z'], self.DataFrame['intensity'])
        plt.xlabel('m/Z')
        plt.ylabel('Intensity (a.u.)')
        plt.show()

## Assigns each folder as a Run, should receive data about how to assign samples
class Run:
    def __init__(self,folder_path, average_samples = False):
        self.folder_path = folder_path
        
        #init should return list of datafiles or list of samples
        #bring up all filenames in folder
        file_list = get_txt_files(self.folder_path)
        ##check what all filenames in folder have in common
        name = os.path.commonprefix(file_list)
        ###assign to Run.name
        if name:
            self.name = name
        else:
            print('No common prefix name for files in this run. Check files in '+self.folder_path)

        ##MOVE THIS TO DATAFILE CLASS
        '''
        #OPTION 1, use maps folder assigned earlier
        maps_path = project_path+maps_folder_name+'/'
        identifier_csv_path = ''
        if os.path.exists(maps_path):
            #if filename matches run name
            maps = get_all_files(maps_path)
            for csv in maps:
                #better method needed, not index based (use split and check if lists are equivalent?), removed extra underscore from end of name
                if self.name[:-1] in csv:
                  identifier_csv_path = maps_path+csv #filename, once found
        #other methods -- notes, default based on Liz's data
                                  
        df = pd.read_csv(identifier_csv_path, header=None)
        #df is assumed to have the first column be cell number, second column be sample name (CAL for calibrant)
        second_column = df.iloc[:,1]
        filenames_to_create = list(set(second_column))
        #remove nan values, ignore CAL
        filenames_to_create = [x for x in filenames_to_create if (str(x)!='nan' and str(x)!='CAL') ]
        

        master_list = []
        #rethink how to structure the data --- TRY TO MAKE DICT IN RUNS AND PASS IN LIST OF RUNS TO SAMPLE
        for suffix in filenames_to_create:
            new_list = []
            for i in range(len(second_column)):
                if second_column[i] == suffix:
                    new_list.append(str(df.iloc[i,0]))
            master_list.append(new_list)
        self.master_dict = dict(zip(filenames_to_create, master_list))
        '''
        if average_samples:
            return [Sample(
                list_of_datafiles=None, output_name=None, folder_path=None, output_path=None
            )]
        return [Datafile(txt_filepath) for txt_filepath in file_list]
            
    def get_data(self):
        #create new directory for sample data
        output_path = self.folder_path+self.name+'Samples/'
        if not os.path.exists(output_path):
            os.makedirs(output_path)
            
        #Associate keys with sample names, make file output
        for key in self.master_dict.keys():
            files_to_use = self.master_dict.get(key)
            Sample(files_to_use, key, self.folder_path, output_path).output_to_file()
            #export sample as txt
                        
                

##Assigns quantities based on MALDI filename convention, more useful if sample names aren't already assigned
class Title:
    def __init__(self, path):
        self.filename = path.removesuffix('.txt')
        meta = self.filename.split(sep='_')
        #assumes the naming convention separates with underscores, putting row and column last
        self.row = meta[-1][0]
        self.column = int(meta[-1][1:])
        self.plate = meta[-1]
        self.date = meta[0]
        self.calibrant = False
        #run 96-well plate structure : Calibrant is every third column starting with the second, and 
        calibrant_rows = ['B', 'E', 'H', 'K', 'N', 'P']
        if (self.row in calibrant_rows) and ((self.column % 3) == 2):
            self.calibrant = True
            #Associate it with other datafiles around the calibrant?
    pass

'''
#ACTUAL PROGRAM - RUN ALL DATA IN THE ASSIGNED PROJECT FOLDER



'''



##Might need for default sample name assignment
'''
#arranging triplicate combinations into samples as per Liz's use
def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in range(ord(c1), ord(c2)+1):
        yield chr(c)

def sample_list(grouping_file_path = None):        

        #default algorithm for full calibrant set, every three columns, plus calibrant
        rows = char_range('A', 'P')
        cols = range(1,25)
        calibrant_rows = ['B', 'E', 'H', 'K', 'N', 'P']
        master_list = []

    
        ###define based on master dict from file
        for row in rows:
            
            for i in range(3):
                sample = []
                for col in cols:
#remove all sets which are calibrant
                    if not ((row in calibrant_rows) and (i%3 == 2)):
                        if (col%3) ==i:
                            sample.append(row+str(col))
            #check if list is not empty
            
                if sample:
                    #only using 1-18 in groups of 3
                    #fix later for full sets
                    master_list.append(sample[0:3])
                    master_list.append(sample[3:6])
        return master_list
    
    #
    

sample_list(grouping_file_path)
'''
'''
#FULL IMPLEMENTATION HERE FOR NOW!!!!********
#LOAD CSV file with known peaks for comparison
df = pd.read_csv(known_peaks_path, header=10)
df_vals = df[df.columns[4:-2]]
df_meta = df[df.columns[[0,1,2,3,-2,-1]]]
print(df_meta)
#indices should match between list and df_meta to correctly assign fish names
list_of_known_peaks = [np.asarray(df_vals.loc[i].dropna()) for i in range(len(df_vals))]

##get all folders in project directory            
folder_list = [ f.path for f in os.scandir(project_path) if f.is_dir() ]
##ignore mapping folder
###NOTE TO SELF: SHOULD CHECK IF IT EXISTS FIRST
folder_list.remove(project_path+maps_folder_name)

##Each folder should represent a single run, with text files for each well in the plate 
for run_folder in folder_list:
    #Performs averages and ouputs to file (see get_data function in Run class)
    for datafile in get_txt_files(run_folder):
        obj = Datafile(run_folder+'/'+datafile)
        file_peaks = obj.get_peak_list()
        print('Run: {}    Cell: {}    Sample: {} \n'.format(obj.run_name, obj.cell, obj.sample_id))
        #list the length of each matches list to align indices
        match_lengths =[]
        for known_peaks in list_of_known_peaks:
            matches, remainder1, remainder2 = list_similarity_operation(file_peaks, known_peaks)
            match_lengths.append(len(matches))
        #assuming this puts indices in reverse order of highest
        indices = np.argpartition(match_lengths, -3)[-3:]
        print('Highest number of peak matches:')
        for index in indices:
            top_val = 3-np.where(indices==index)[0]
            top_spec = df_meta.loc[index]['Species']
            top_len = match_lengths[index]
            print('#{}: {}, {} matches'.format(top_val, top_spec, top_len))
        print('Datafile analysis complete. \n')
        
    
##If completed with no errors, acknowledge    
print('Project averaging completed from '+project_path)

'''
#PROGRAM OUTLINE - PERSONAL USE

##access folder as a run
##each run is ordered into samples
##each sample is averaged from n datafiles