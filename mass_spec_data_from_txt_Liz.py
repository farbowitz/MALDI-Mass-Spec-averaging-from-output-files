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
from scipy.ndimage import gaussian_filter1d

#ATTEMPTING LOGGING
import logging
import sys

from matplotlib.ticker import FormatStrFormatter
from Datafile import Datafile

logger = logging.getLogger(__name__)
stream_handler = logging.StreamHandler(sys.stdout)
logger.addHandler(stream_handler)

# FUNDAMENTAL VARIABLE INPUTS
known_peaks_path="C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Mass_spec_peaks_fish.csv"
project_path = "C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Usable Spectra July 2022/"
## OPTIONAL- SAMPLE NAMING METHOD: assign a folder to keep CSV files (titled with the run name) which map column 1 wells (e.g. A1, B3, etc.) to sample names (e.g. sample 2, coppergate_fish). Use CAL for calibrant wells. 
maps_folder_name = 'maps'
## default - adds folder for sample data under projet directory

###TEST DATA
test_folder = project_path+'20220711_LMQ_Test Text Files/'
test_file = '20220711_LMQ_Test_300622_A3.txt'
test_path = test_folder+test_file


#USEFUL FUNCTIONS

def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in range(ord(c1), ord(c2)+1):
        yield chr(c)


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

def get_map_pd(folder=test_folder):
    map_folder = project_path+maps_folder_name+'/'
    if os.path.exists(map_folder):
        #if filename matches run name
        maps = get_all_files(map_folder)
        for map in maps:
            #better method needed, not index based (use split and check if lists are equivalent?), removed extra underscore from end of name
            if 'LMQ_Test' in map:
                map_path = map_folder+map #filename, once found
        
            
        logger.warning('No map found containing {} in {}.'.format(maps_folder_name, 'LMQ_test'))
            
    else:
        logger.warning('No maps folder found. Creating one at {}'.format(map_folder))
        #create folder
        #add default.csv file to it
        map_path=None
    df = pd.read_csv(map_path, names=['cell id', 'sample id'], header=None)
    #df is assumed to have the first column be cell number, second column be sample name (CAL for calibrant)
    df = df.dropna()
    #remove nan values, ignore 

    return df

#not working
def batch_convert(folder=test_folder):
    file_list = get_txt_files(folder)
    df = get_map_pd()

    #create Renamed_files in folder
    directory = test_folder+'Renamed_Files/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    for file in file_list:
        name = file.replace('.txt', '')
        cell_name = name.split(sep='_')[-1]

        sample_name = df.loc[df['cell id']==cell_name, 'sample id'].item()
        new_name = name.replace(cell_name, sample_name)
        chars = ['A', 'B', 'C', 'D', 'E', 'F']
        for index in range(len(chars)):
            try_name = new_name + '_' + chars[index] + '.txt'
            if not os.path.exists(directory+new_name):
                new_name = try_name
                break
            print(folder+file)
            print(directory+new_name)
            with open(folder+file,'r') as firstfile, open(directory+new_name,'w') as secondfile:
                # read content from first file
                for line in firstfile:
                    print(line)
                    # append content to second file
                    secondfile.write(line)
                secondfile.close()
                firstfile.close()

#LOAD CSV file with known peaks for comparison
def load_known_peaks(self, known_peaks_path="C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Mass_spec_peaks_fish.csv"):
    df = pd.read_csv(known_peaks_path, header=10)
    df_vals = df[df.columns[4:-2]]
    df_meta = df[df.columns[[0,1,2,3,-2,-1]]]

    #indices should match between list and df_meta to correctly assign fish names
    list_of_known_peaks = [np.asarray(df_vals.loc[i].dropna()) for i in range(len(df_vals))]
    return list_of_known_peaks, df_meta







# CLASSES




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

def list_similarity_operation(list1, list2, allowable_distance=1):
    match_list = []
    non_matches1 = []
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

#ATTEMPT ONE: LIZ THINKS I CAN JUST DO ANYTHING

#test_run = Run(folder_path = project_path+'30062022_LMQ_Test Text Files/')

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

#FULL IMPLEMENTATION HERE FOR NOW!!!!********
#keep as outside function or class function?
def print_species_matches_from_datafile(path):
        obj = Datafile(path)
        file_peaks = obj.get_peak_list()
        print('Sample: {}'.format(obj.sample_id))
        list_of_matches, match_lengths, df_meta = obj.list_matches_to_fish_data(file_peaks)
        #list the length of each matches list to align indices

        #assuming this puts indices in reverse order of highest
        indices = np.argpartition(match_lengths, -3)[-3:]
        print('Highest number of peak matches:')
        for index in indices:
            top_val = 3-np.where(indices==index)[0]
            top_spec = df_meta.loc[index]['Species']
            top_len = match_lengths[index]
            print('#{}: {}, {} matches @ {}'.format(top_val, top_spec, top_len, list_of_matches[index]))
        print('\n')
        logger.info('Datafile analysis complete. \n')


def comparison_metric(x1_peaks, y1_peaks, x2_peaks, y2_peaks):
    
    '''
    minimum = 0.1
    #adjust peak list 
    x1_peaks_major = x1_peaks[np.where(y1_peaks > minimum)]
    y1_peaks_major = y1_peaks[np.where(y1_peaks > minimum)]
    x2_peaks_major = x2_peaks[np.where(y2_peaks > minimum)]
    y2_peaks_major = y2_peaks[np.where(y2_peaks > minimum)]
    '''
    
    matched_peaks = []
    for peak in x1_peaks:
        index1 = int(np.where(x1_peaks == peak)[0])
        index2 = int(find_nearest_index(x2_peaks, peak))
        x1 = x1_peaks[index1]
        x2 = x2_peaks[index2]
        y1 = y1_peaks[index1]
        y2 = y2_peaks[index2]
        #print('y2: {} y1: {}   x2: {}  x1:{}'.format(y2,y1,x2,x1))
        #print('y_distance: {}    x_distance: {}'.format(np.abs(y2-y1), np.abs(x2-x1)))
        if np.abs(x2-x1) < 0.1:
            matched_peaks.append((x1, x2, y1, y2))
        #euclidean_distance = np.sqrt((x2-x1)**2+(y2-y1)**2)
        #np.append(distances, euclidean_distance)
    
    #return distances
    return matched_peaks

###Used to limit plots to grouping the closest to square
def squarest_dims(number):
    #find nearest square larger than number
    root = np.ceil(np.sqrt(number))
    #can it fit into n by (n-1) dimensions?
    if number < ((root-1)*root):
        return (root-1), root
    else:
        return root, root
    
def plot_within_subplots(obj, df, axes, index_x, index_y, x1, x2, y1, y2):
    axes[index_x,index_y].plot(df['m/Z'], df['intensity']/obj.DataFrame['intensity'].max())
    axes[index_x,index_y].annotate(round(x1,2), xy=(x1, y1))
    axes[index_x,index_y].annotate(round(x2,2), xy=(x2, y2))
    axes[index_x, index_y].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    #axes[index_x,index_y].title(np.average([x1,x2]))

def full_comparison(obj1, obj2):
    name1 = obj1.sample_id
    name2 = obj2.sample_id
    x1_peaks, y1_peaks = obj1.peak_coords()
    x2_peaks, y2_peaks = obj2.peak_coords()
    print('\n Comaprison between {} and {}:'.format(name1,name2))
    obj1.plot_peaks()
    obj2.plot_peaks()
    matched_peaks = comparison_metric(x1_peaks, y1_peaks, x2_peaks, y2_peaks)
    print('\n Number of matches: {}'.format(len(matched_peaks)))
    n_of_graphs = len(matched_peaks)
    nrows, ncols = squarest_dims(n_of_graphs)
            
    fig, axes = plt.subplots(nrows=int(nrows), ncols=int(ncols))
    for (x1,x2, y1, y2) in matched_peaks:
        avg = np.average([x1, x2])
        dist = 0.5
        frame_max = avg+dist
        frame_min = avg-dist
        df1 = obj1.clip_data(obj1.DataFrame, frame_min, frame_max)
        df2 = obj2.clip_data(obj2.DataFrame, frame_min, frame_max)
        
        idx = matched_peaks.index((x1,x2,y1,y2))
        index_x = int(np.floor(idx/ncols))
        index_y = int(idx % ncols)
        plot_within_subplots(obj1, df1, axes, index_x, index_y, x1, x2, y1, y2)
        plot_within_subplots(obj2, df2, axes, index_x, index_y, x1, x2, y1, y2)
        

    fig.suptitle('{} and {} peak matches'.format(name1,name2))
    fig.supxlabel('m/Z')
    fig.supylabel('Intensity (a.u.)')
    plt.tight_layout()
        
    plt.show()

#test code for datafile class
test_path = 'C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Usable Spectra July 2022/LMQ116_A.txt'
obj1 = Datafile(test_path)
test_path=test_path.replace('LMQ116_A', 'LMQ052_B')
obj2 = Datafile(test_path)

obj1.plot_peaks()
obj2.plot_peaks()
#x1_peaks, y1_peaks = obj1.peak_coords()
#x2_peaks, y2_peaks = obj2.peak_coords()
#print(comparison_metric(x1_peaks, y1_peaks, x2_peaks, y2_peaks))
#full_comparison(obj1, obj2)


    

'''
##get all folders in project directory            
folder_list = [ f.path for f in os.scandir(project_path) if f.is_dir() ]
##ignore mapping folder
###NOTE TO SELF: SHOULD CHECK IF IT EXISTS FIRST
folder_list.remove(project_path+maps_folder_name)


##Each folder should represent a single run, with text files for each well in the plate 
for run_folder in folder_list:
'''  
'''
#Performs averages and ouputs to file (see get_data function in Run class)
for datafile1 in get_txt_files(project_path):
    for datafile2 in get_txt_files(project_path):
        obj1 = Datafile(project_path+'/'+datafile1)
        obj2 = Datafile(project_path+'/'+datafile2)

        full_comparison(obj1, obj2)
        
'''      
       
        
    
##If completed with no errors, acknowledge    
print('Project averaging completed from '+project_path)


#PROGRAM OUTLINE - PERSONAL USE

##access folder as a run
##each run is ordered into samples
##each sample is averaged from n datafiles