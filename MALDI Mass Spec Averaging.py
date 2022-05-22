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
from scipy.signal import find_peaks
from Mass_Spec_Compare import peak_list
from lmfit.models import LorentzianModel, Model, GaussianModel, SplitLorentzianModel, MoffatModel, VoigtModel, PseudoVoigtModel, Pearson7Model, StudentsTModel, LinearModel

# FUNDAMENTAL VARIABLE INPUTS

project_path = "C:/Users/Daniel/Desktop/Liz_MALDI_runs/"
## OPTIONAL- SAMPLE NAMING METHOD: assign a folder to keep CSV files (titled with the run name) which map column 1 wells (e.g. A1, B3, etc.) to sample names (e.g. sample 2, coppergate_fish). Use CAL for calibrant wells. 
maps_folder_name = 'maps'
## default - adds folder for sample data under projet directory
default_csv_url = 'https://github.com/farbowitz/MALDI-Mass-Spec-averaging-from-output-files/blob/main/Default%20-%20Sheet1.csv'
##imports database of peaks as list of lists, see Mass_Spec_Compare
peak_list = peak_list()
##define region of interest for m/Z values
region_of_interest = [850, 3000]


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

##Might need for default sample name assignment

##Finds the index of a list based on closest value to given number
def find_nearest_dex(array, number, direction=None): 
    idx = -1
    if direction is None:
        ser = np.abs(array-number)
        idx = ser.get_loc(ser.min())
    elif direction == 'backward':
        _delta = number - array
        _delta_positive = _delta[_delta > 0]
        if not _delta_positive.empty:
            idx = _delta.get_loc((_delta_positive.min()))
    elif direction == 'forward':
        _delta = array - number
        _delta_positive = _delta[_delta >= 0]
        if not _delta_positive.empty:
            idx = _delta.get_loc(_delta_positive.min())
    return idx

## Graphing Functions
###get points around center
def interval(dataX, dataY, center, deviation):
    bottom = center-deviation
    top = center+deviation
    newX = dataX[bottom:top]
    newY = dataY[bottom:top]
    return newX, newY

###plot using models given via lmfit.models library, e.g. GaussianModel()
def plot_with_model(X, Y, model, mname):
    params = model.guess(Y, x=X, cen=X[8], amp=max(Y)-min(Y))
    result = model.fit(Y, params, x=X)
    aic = result.aic
    plt.scatter(X, Y, c='black')
    plt.plot(X, result.best_fit, '--', label=mname+' AIC: '+str(np.around(aic,3)))
    plt.xlabel('m/Z')
    plt.ylabel('intensity(a.u.)')
    #leaves out plt.plot(), intended to be looped prior to plt.plot()

### break data into smaller sections by values rather than index
def subset_by_x_values(dataX, dataY, bottom, top):
    bottom_index = find_nearest_dex(dataX, bottom)
    top_index = find_nearest_dex(dataX, top)
    newX = dataX[bottom_index:top_index]
    newY = dataY[bottom_index:top_index]
    return newX, newY

### define number of intervals to break region into





# CLASSES

## individual txt files from each run are a datafile
class Datafile:
    def __init__(self, txt_filepath):
        #import data into pandas based on filename referenced
        df = pd.read_csv(txt_filepath, sep="\t", header=0)
        df.columns = ['m/Z', 'intensity']
        #change m/Z column to be index
        df.set_index('m/Z', inplace=True)
        self.DataFrame = df
        #return none if no data
            
    def identify_peaks(self):
        return None
    def clean_noise(self):
        return None
    #internal function to plot data
    def plot(self):
        plt.plot(self.DataFrame['m/Z'], self.DataFrame['intensity'])
        plt.xlabel('m/Z')
        plt.ylabel('Intensity (a.u.)')
        
        

## Creates a Sample from one or multiple datafiles
class Sample:
    def __init__(self, list_of_datafiles, output_name, folder_path, output_path):
        self.size = len(list_of_datafiles)
        self.list = list_of_datafiles
        self.name = output_name
        self.folder_path = folder_path
        self.output_path = output_path
        #make empty file of column length, *SHOULD ASSIGN USING ARRAY LENGTH, UNLESS ALL MALDI RUNS ARE SAME SIZE*, len(df) in pandas 
        run_file_list = get_txt_files(self.folder_path)
        #self.unmodified = True
        #default maldi file length should be 107875
        #running_total = [0] * 107875
        list_as_string = '-'.join(self.list)
        '''Commented to go easier on computer
        #check if files making up the sample are same size
        array_sizes = []
        for file in run_file_list:
            filename = Title(file)
            if (str(filename.row)+str(filename.column)) in self.list:
                array_sizes.append(len(np.asarray(Datafile(self.folder_path+file).DataFrame['intensity'])))
        if all([x==array_sizes[0] for x in array_sizes]) and array_sizes:
            #adjust running_total list to be length of sample if not 
            running_total = [0] * int(array_sizes[0])
        else:
            print('Array sizes not equal, list not passed')'''
            
        #TRY USING PD DATAFRAMES IN DICT {}
        sample_dataframes = {}
        for file in run_file_list:
            #check if cell is in list_of_datafiles, SHOULD APPEND TO DICT
            filename = Title(file)
            cell_name = str(filename.row)+str(filename.column)
            if cell_name in self.list:
                sample_dataframes[cell_name] = Datafile(self.folder_path+file).DataFrame
        #check if any dfs exist in sample_dataframes
        if not sample_dataframes:
            print('No files found containing '+list_as_string)
            return None
        #use first dataframe as reference for dataframe shape and index
        df1_name = list(sample_dataframes.keys())[0]
        df1 = sample_dataframes[df1_name]
        reference_index = df1.index
        reference_shape = df1.shape
        #compare all dataframes, check for inconsistencies
        for cell_name in sample_dataframes:
            dfx = sample_dataframes[cell_name]
            #check index is the same
            ##PROBLEM: MOST INDICES DON'T MATCH, THIS POSES AN ISSUE FOR AVERAGING
            if not (reference_index.equals(dfx.index)):
                print('Indices don\'t match at '+str(cell_name))
                return None
            ###SOLUTION 1: PUT ALL DATA TOGETHER AND ROLL WITH IT
            ###SOLUTION 2:BIN DATA AND AVERAGE
            #check for same dimensions
            if not (reference_shape == dfx.shape):
                print('Dataframe sizes don\'t match. '+df1_name+' gives '+df1.shape+' and '+str(cell_name)+' gives '+ dfx.shape)
                return None
            
        #check that files were found
        ##NOTE TO SELF: STOP SWITCHING BETWEEN LIST AND ARRAYS SO MUCH -- FIGURE OUT HOW TO BE CONSISTENT, USE PANDAS DFS?
        print('Averaging ' + self.name)
        df = pd.Panel(sample_dataframes).mean(axis=0)
        print(df)
        self.DataFrame = df
        
    def output_to_file(self):
        if not self.DataFrame:
            return None
        else:
            self.DataFrame.to_csv(r''+self.output_path+self.name+'.txt', sep='\t', index=True, header=False)
        
    def identify_peaks(self):
        return None
    
    def clean_noise(self):
        return None        
    #Returns keyError -- data type issue?        
    def plot(self):
        plt.plot(self.DataFrame['m/Z'], self.DataFrame['intensity'])
        plt.xlabel('m/Z')
        plt.ylabel('Intensity (a.u.)')
    
    #rethink these functions, avoid replicating methods
    def find_peaks(self):
        X = self.DataFrame.index
        y = self.DataFrame.intensity
        #prominence is trough to peak tolerance for peak identification, set to a tenth of data's range
        prominence = (np.max(y)-np.min(y))/10
        #others --- height, width, tolerance
        peaks = find_peaks(y, prominence=prominence)
        peak_index = peaks[0]

        plt.plot()
        y_peaks = [y[index] for index in peak_index]
        x_peaks = [X[index] for index in peak_index]
        plt.scatter(x_peaks, y_peaks, c='orange')
        plt.show()
        self.x_peaks = x_peaks
        self.y_peaks = y_peaks
        
    def infer_peak(self):
        X = self.DataFrame.index
        Y = self.DataFrame.intensity
        fish_name_list = []
        for item in peak_list:
            fish_name_list.append(item[1])
        #severly confusing data types and methods here
        fish_names = ['Atlantic salmon']
        indexer = [(any(fish_names) == name) for name in fish_name_list]
        salmon_peaks = np.asarray(peak_list)[np.asarray(indexer)]
        for i in range(len(salmon_peaks)):
            source_name = salmon_peaks[i][2]
            source_peaks = salmon_peaks[i][0]
            ''' # Looks at each individual peak
            for peak in source_peaks:
                Xn, Yn = interval(X,Y, find_nearest_dex(X, peak), 100)
                plt.plot(Xn,Yn)
                plt.axvline(x=peak, c='orange', label='m/Z = '+str(peak)+', source(s): '+source_name)
                plt.legend(loc='best')
                plt.title('Sample '+self.name+ ' Local Peaks')
                plt.show()'''
                
            #check between 850 and 3000, signal-to-noise ratio = 6
            
            plt.plot(X,Y)
            plt.vlines(source_peaks, 0, 5000, colors='orange')
            plt.show()
        
        #should attempt looking for half max widths (either in terms of m/Z or by # of data points), check for consistency
        
    def peak_characterization(self):
        #data
        X = self.DataFrame.index
        Y = self.DataFrame.intensity
        '''#are step sizes all the same?
        delta_xs = [(X[i]-X[i-1]) for i in range(1,len(X))]
        print(delta_xs[0])
        print(np.all(np.asarray(delta_xs)) == delta_xs[0])
        #false, so what does delta_xs look like
        print(delta_xs[57000:57100])
        #not the same, but close! check:
        print(all([(0.017<dx<0.042) for dx in delta_xs]))
        #delta_xs appears to increase from ~0.018 to 0.041
        '''
        #playing around looking for individual peaks
        #major prominence around 1220.6, index= 20521
        
        #sample_peak = 20521
        #num_points around peak: ~5 to half max
        #delta_s = 8
        
        sample_X, sample_Y = interval(X, Y, 20521, 8)
        #recenter around maximum
        max_index = np.argmax(sample_Y)
        sample_X, sample_Y = interval(X,Y,20521+max_index-7, 8)
        #reapproriated from my phy6011 lab 3 collab file
        models = [PseudoVoigtModel(), GaussianModel(), LorentzianModel()]
        mnames = ['PseudoVoigt', 'Gaussian', 'Lorentzian']
        for i in range(len(models)):
            plot_with_model(sample_X, sample_Y, models[i], mnames[i])
        plt.legend(loc='best')
        plt.show()
        
        
        
## Assigns each folder as a Run, should receive data about how to assign samples
class Run:
    def __init__(self,folder_path):
        self.folder_path = folder_path
        
        #bring up all filenames in folder
        file_list = get_txt_files(self.folder_path)
        ##check what all filenames in folder have in common
        name = os.path.commonprefix(file_list)
        #take away parts not date (first) and cell number (last)
        name = '_'.join(name.split(sep='_')[1:-1])
        ###assign to Run.name
        if name:
            self.name = name
        else:
            print('No common prefix name for files in this run. Check files in '+self.folder_path)

        #OPTION 1, use maps folder assigned earlier
        maps_path = project_path+maps_folder_name+'/'
        identifier_csv_path = maps_path+'Default.csv'
        if os.path.exists(maps_path):
            #if filename matches run name
            maps = get_all_files(maps_path)
            for csv in maps:
                #better method needed, not index based (use split and check if lists are equivalent?), removed extra underscore from end of name
                if self.name in csv:
                  identifier_csv_path = maps_path+csv #filename, once found
                #otherwise, use default csv map from github url
        print('Using sample names from '+identifier_csv_path)                          
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
            
    def get_data(self):
        #create new directory for sample data
        output_path = self.folder_path+self.name+'Samples/'
        if not os.path.exists(output_path):
            os.makedirs(output_path)
            
        #Associate keys with sample names, make file output
        for key in self.master_dict.keys():
            files_to_use = self.master_dict.get(key)
            #export sample as txt
            #find peaks and plot
            Sample(files_to_use, key, self.folder_path, output_path)
    
                        
                

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


#ACTUAL PROGRAM - RUN ALL DATA IN THE ASSIGNED PROJECT FOLDER

##get all folders in project directory            
folder_list = [ f.path for f in os.scandir(project_path) if f.is_dir() ]
##ignore mapping folder
###NOTE TO SELF: SHOULD CHECK IF IT EXISTS FIRST
folder_list.remove(project_path+maps_folder_name)

##Each folder should represent a single run, with text files for each well in the plate 
for run_folder in folder_list:
    #Performs averages and ouputs to file (see get_data function in Run class)
    Run(run_folder+'/').get_data()
    
##If completed with no errors, acknowledge    
print('Project averaging completed from '+project_path)




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


#PROGRAM OUTLINE - PERSONAL USE

##access folder as a run
##each run is ordered into samples
##each sample is averaged from n datafiles

#Possible Additions
##Directly convert from XML output files (no R necessary!)
##Statistical test prior to averaging -- check samples are same (Check if Kolmogorov-Smirnov Test is appropriate one)
##Peak detection
##Peak comparisons with known data
###Statistical test for bone identification and degree of certainty
##Plot data