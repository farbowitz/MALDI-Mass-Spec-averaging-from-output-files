# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 12:23:58 2022

@author: Daniel

Machine Learning Models For Zoom

1st idea: Use list of known peaks for species, make pd DataFrame with all peaks listed, 0 for not present, 1 for present

Information theory?
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mass_spec_data_from_txt_Liz import Datafile

import logging
import sys

l = logging.getLogger(__name__)
stream_handler = logging.StreamHandler(sys.stdout)
l.addHandler(stream_handler)


known_peaks_path="C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Mass_spec_peaks_fish.csv"

class LearningModel:
    def __init__(self, key=known_peaks_path):
        self.key = key
        #self.peaks, self.identifier = self.import_data()
        self.peak_values = self.get_unique_peaks()
        self.candidate_species = self.import_data().columns
        print(self.candidate_species)
        self.matrix = self.create_new_dataframe()
        self.dataset, self.labels = self.initial_train_test_from_dataframe()
        self.seed = np.random.randint(50)

        #look up / narrowing for genus traits?
        
        
        pass
    def import_data(self):
        '''
        

        Returns
        -------
        df_vals : pd.DataFrame of peaks with species names as columns
            DESCRIPTION.

        '''
        key = self.key
        #df_meta = df[df.columns[[0,1,2,3,-2,-1]]]
        #in this case, casts data into species being column names
        df = pd.read_csv(key, header=10)
        df_vals = df[df.columns[3:-2]]
        df_vals = df_vals.set_index(['Species'])
        df_vals.index = df_vals.index.where(~df_vals.index.duplicated(), df_vals.index + '_dp')
        #separate data with same species name
        
        df_vals = df_vals.T
        return df_vals
    
    def get_unique_peaks(self):
        df = self.import_data()
        array_nans_and_unsorted = pd.unique(df[df.columns].values.ravel())
        array_unsorted = array_nans_and_unsorted[~np.isnan(array_nans_and_unsorted)]
        array = sorted(array_unsorted)
        return array
    
    def create_new_dataframe(self):
        #RESOURCE INTENSIVE, must be better implementation
        df = self.import_data()
        species = df.columns
        peak_values = self.get_unique_peaks()
        analysis_df = pd.DataFrame(index=species, columns=peak_values)
        for specie in species:
            for peak in peak_values:
                if peak in (pd.unique(df[specie].values.ravel())):
                    analysis_df.loc[specie, peak] = 1
                else:
                    analysis_df.loc[specie, peak] = 0
        return analysis_df
    
    def initial_train_test_from_dataframe(self):
        df = self.matrix
        init_dataset = []
        init_labels = []
        for specie in df.index:
            init_dataset.append(list(df.loc[specie]))
            init_labels.append(specie)
        return init_dataset, init_labels

    def amplified_initial_dataset(self, n=50):
        init_dataset, init_labels = self.initial_train_test_from_dataframe()
        return n*init_dataset, n*init_labels
    
    def import_test_data(self):
        '''
        Imports peak info from runs with known species.

        Returns
        -------
        None.

        '''
        pass

    def import_single_run(self, file):
        run_peaks = Datafile(file).get_peak_list()
        species_id = Datafile(file).get_species_id()
        return run_peaks, species_id

    def characterize_run_to_list(self, file):
        run_peaks, species_id = self.import_single_run(file)
        if species_id is not None:
            l.info('Species ID detected, adding to list.')
            characterization_list = []
            for peak in self.peak_values:
                sample_array = np.asarray(run_peaks)
                peak_distances = np.abs(sample_array-peak)
                if np.any(peak_distances<0.5):
                    characterization_list.append(1)
                else:
                    characterization_list.append(0)
            #double check size is correct
            if len(characterization_list) == len(self.peak_values):
                #update model's dataset and labels
                self.dataset.append(characterization_list)
                self.labels.append(species_id)
            else:
                l.warning('Size difference detected.')
        else:
            l.info('Species not IDed yet.')

    def characterize_unIDed_sample(self, file):
        run_peaks, species_id = self.import_single_run(file)
        characterization_list = []
        for peak in self.peak_values:
            sample_array = np.asarray(run_peaks)
            peak_distances = np.abs(sample_array-peak)
            if np.any(peak_distances<0.5):
                characterization_list.append(1)
            else:
                characterization_list.append(0)
        if len(characterization_list) == len(self.peak_values):
            #update model's dataset and labels
            return characterization_list
        else:
            l.warning('Size difference detected.')
            return None
    
    def classify_unIDed_sample(self, file):
        c_list = self.characterize_unIDed_sample(file)
        train_data, test_data, train_labels, test_labels = self.get_vals()
        classifier = self.decision_tree().fit(train_data, train_labels)
        output = classifier.predict(np.asarray(c_list).reshape(1,-1))
        return output

    ##PREPPING MODEL TRAINING/TEST DATA

    def get_vals(self):
        from sklearn.model_selection import train_test_split
        return train_test_split(self.dataset, self.labels, random_state=self.seed, test_size=0.2)
    
    def get_test_vals(self):
        _, test_data, _, test_labels = self.get_vals()
        return test_data, test_labels

    def get_train_vals(self):
        train_data, _, train_labels, _ = self.get_vals()
        return train_data, train_labels

    ##MODELS

    def decision_tree(self, i=75):
        from sklearn.tree import DecisionTreeClassifier
        classifier = DecisionTreeClassifier(random_state=self.seed, max_depth=i)
        return classifier

        pass
    
    def forest(self):
        from sklearn.ensemble import RandomForestClassifier
        classifier = RandomForestClassifier(n_estimators = len(self.dataset))
        return classifier
    
    def naive_bayes(self):
        #from sklearn.naive_bayes import *
        pass
    def svm(self):
        #from sklearn.svm import *
        pass
    
    
    ##RETURNS BASED ON MODELS
    
    def classifier_to_score(self, classifier):
        train_data, test_data, train_labels, test_labels = self.get_vals()
        classifier.fit(train_data, train_labels)
        return classifier.score(test_data, test_labels)

    def plot_scores(self, max_value=100):
        scores = []
        values = []
        for i in range(len(self.candidate_species)):
            scores.append(self.classifier_to_score(self.decision_tree(i=i+1)))
            values.append(i+1)
        plt.plot(values,scores)
        plt.title('Tree classifier scores (out of 1) based on tree depth')
        plt.show()


        #indices should match between list and df_meta to correctly assign fish names
        #list_of_known_peaks = [np.asarray(df_vals.loc[i].dropna()) for i in range(len(df_vals))]
        #return list_of_known_peaks, df_meta

test_path = 'C:/Users/Daniel/Desktop/Programming/Mass Spec Project/Usable Spectra July 2022/C08_A.txt'
model = LearningModel()
id = model.classify_unIDed_sample(test_path)
print(id)
a = np.asarray(model.characterize_unIDed_sample(test_path))
b = np.asarray(model.matrix.loc[id])
print(len(a))
print(np.sum(a==b))