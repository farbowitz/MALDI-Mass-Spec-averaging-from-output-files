# -*- coding: utf-8 -*-
"""
Created on Sun May 22 00:29:09 2022

@author: Daniel
"""

'''
Testing basics for MALDI Mass Spec Averaging
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sample_path = 'C:/Users/Daniel/Desktop/Liz_MALDI_runs/20220429_Coppergate_1 Text Files/20220429_Coppergate_1_'
sample_names = ['A2', 'A5', 'A8']
sample_dfs = []
for name in sample_names:
    df = pd.read_csv(sample_path+name+'.txt', sep="\t", header=0)
    df.columns = ['m/Z', 'intensity']
    df.set_index('m/Z', inplace=True)
    sample_dfs.append(df)
    
'''
#from dict earlier
df1_name = list(sample_dfs.keys())[0]
df1 = sample_dfs[df1_name]
reference_index = df1.index
reference_shape = df1.shape
df2_name = list(sample_dfs.keys())[1]
df2 = sample_dfs[df2_name]
second_index = df2.index
'''
###SOLUTION 1: PUT ALL DATA TOGETHER AND ROLL WITH IT
'''df_combined = pd.concat([df1,df2])
plt.scatter(df1.index, df1.intensity)
plt.scatter(df2.index, df2.intensity)
#WTF scale isn't even remotely the same -- A1 max is 14,000, A7 max is ~4,000, peaks appear the same though
#Makes combining data impossible
plt.show()
plt.scatter(df_combined.index, df_combined.intensity)
plt.show()
'''
###SOLUTION 2:BIN DATA AND AVERAGE
#print(df1.index.min())
#print(df2.index.min())
##Noticed deltas are monotonically increassing, so
'''delta_min = df1.index[1]-df1.index[0]


print(delta_min)
'''


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
'''
top_index = find_nearest_dex(df1.index, 3200)

delta_max = df1.index[top_index]-df1.index[top_index-1]
print(delta_max)


    #deltas slowly increase with index


'''
#try using function to relate, how close can we get?
def find_nearest(array, number, direction): 
    idx = -1
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

'''
natural_deltas = []

adjusted_deltas = []
for i in range(len(reference_index)):
    natural_deltas.append(reference_index[i]-second_index[i])
    adjusted_deltas.append(find_nearest(np.asarray(second_index), reference_index[i], None))
#interesting, may be a limit on this. How close and how far do we get?
print('Deltas from index-based approach')
print('Max:' + str(max(np.abs(natural_deltas))))
print('Min:' + str(min(np.abs(natural_deltas))))
print('Modified Deltas from find_nearest function')
print('Max:' + str(max(adjusted_deltas)))
print('Min:' + str(min(adjusted_deltas)))
'''

#SOLUTION 2 (AGAIN)
def top_step_size(index, top_cutoff_value):
    dex = find_nearest_dex(index, top_cutoff_value, 'forward')
    delta = index[dex]-index[dex-1]
    return delta

def bottom(index, bottom_cutoff_value):
    if bottom_cutoff_value > index.min():
        return bottom_cutoff_value
    else:
        return index.min()

def df_averaging_and_cropping_with_unequal_indices(list_of_dfs, bottom_cutoff_value, top_cutoff_value, rolling_val):
    df_full = pd.concat(list_of_dfs)
    df_sorted = df_full.sort_index()
    df_in_interval = df_sorted[df_sorted.index.to_series().between(bottom(df_sorted.index,bottom_cutoff_value), top_cutoff_value)]
    return df_in_interval.rolling(rolling_val).mean()
df = df_averaging_and_cropping_with_unequal_indices(sample_dfs, 1443.2, 1444.2, 30)
plt.scatter(df.index, df.intensity)
plt.show()

'''
#What it looks like for each individual datafile
for i in range(len(sample_dfs)):
    
    sample_df = sample_dfs[i][sample_dfs[i].index.to_series().between(bottom(sample_dfs[i].index,850), 3000)]
    plt.scatter(sample_df.index, sample_df.intensity, label=sample_names[i])
plt.legend()
plt.xlabel('m/Z')
plt.ylabel('Intensity')
plt.show()     
'''

'''
#what do numeric m/Z peaks look like in terms of index? largest peak ~ 1443-1444
print(find_nearest_dex(sample_dfs[0].index, 1443.7))

for i in range(len(sample_dfs)):
    
    sample_df = sample_dfs[i].reset_index()
    sample_df = sample_df[sample_df.index.to_series().between(29900, 30030)]
    plt.scatter(sample_df.index, sample_df.intensity, label=sample_names[i])
plt.legend()
plt.xlabel('Index Number')
plt.ylabel('Intensity')
plt.show()     
'''

###Rolling method appears to be okay. Couple things to check first:
    #same number of datapoints
    #1) is domain (x-range) for all cells the same?
'''
for i in range(len(sample_dfs)):
    sample_df = sample_dfs[i]
    index = sample_df.index
    name = sample_names[i]
    print('Dataset: {name} '.format(name=name))
    print('Max: {} , Min: {} , Range: {}'.format(max(index), min(index), max(index)-min(index)))'''
    ##Answer: NO, range differs by as much as .8
    #2) do datapoints come from successive datasets? Is it random?
'''sample_dfs2 = []
for i in range(len(sample_dfs)):
    sample_df = sample_dfs[i]
    sample_df['Dataset'] = sample_names[i]
    sample_dfs2.append(sample_df)
arranged_df = pd.concat(sample_dfs2).sort_index()['Dataset']
    ##Answer: Appear to occur successively, e.g. A1, A7, A4... starting a few points in, but then how are ranges different and # of points the same?
    ##Might stick with rolling: graphically looks the same, but I should make sure that assumption is justified somehow

name_seq = ['A4', 'A1', 'A7']
#brute force, sorry computer!
follows_seq = []
for i in range(102,22000):
    follows_seq.append(name_seq[i%3]==arranged_df.iloc[i])
print('Result is {}'.format(all(follows_seq)))
print(arranged_df.index[22000])
#works for 800-940, then what?
print(arranged_df.iloc[24000:24003])
#arrangement changed slightly?
name_seq_2 = ['A7', 'A1', 'A4']
follows_seq_2 = []
for i in range(24000, 27000):
    
    follows_seq_2.append(name_seq_2[i%3]==arranged_df.iloc[i])
print('Result is {}'.format(all(follows_seq_2)))
print(arranged_df.iloc[64000:64010])
#Mostly follows sequence, when does it violate this? Is it small enough rolling(3) is reasonable?
##going to move on with rolling(3) and HOPEFULLY return to this to confirm

    
    
    
    
    #3) what are detla values for each? are they the same?
'''

#Let's try curve analysis using individual datasets rather than the incorrect averages earlier

for i in range(len(sample_dfs)):
    
    sample_df = sample_dfs[i][sample_dfs[i].index.to_series().between(bottom(sample_dfs[i].index,1443.2), 1444.2)]
    plt.scatter(sample_df.index, sample_df.intensity, label=sample_names[i])
plt.legend()
plt.xlabel('m/Z')
plt.ylabel('Intensity')
plt.show()     


    
print('Completed without error.')
