import numpy as np 
import pandas as pd
from sklearn import linear_model, metrics, ensemble, model_selection, preprocessing
from DataLoaderRF import DataLoader
from random import shuffle
from FeatureSelection import *
from utility import *
import matplotlib.pyplot as plt 

files = ['nucleus.fasta.txt', 'cyto.fasta.txt', 'mito.fasta.txt', 'secreted.fasta.txt']
labels = ['nucleus', 'cyto', 'mito', 'secreted']

class GradientBooster(object):
    def __init__(self):
        self.load_data()
        self.data = self.dataLoad.data
        
    def load_data(self):
        self.dataLoad = DataLoader(path='./Data/')
        self.dataLoad.load_files(files, labels)
        self.dataLoad.convert_to_df(verbose=False)
        
    def add_features(self):
        sequence_length(self.data)
        numerate_sequence_and_attributes(self.data)
        expand_attributes(self.data)
        amino_acid_count(self.data)
        amino_acid_percentage(self.data)
        aromatic_percentage(self.data)
        self.data = explode_feature(self.data, 'Amino Acid Percentage')
        self.data = explode_feature(self.data, 'Amino Acid Percentage Start')
        self.data = explode_feature(self.data, 'Amino Acid Percentage End')
        
    def get_labels_and_data(self):
        features = ['Sequence Length Norm', 'Sequence Weights_mean', 'Sequence Weights_std', 'Sequence Isoelectric_mean',
            'Sequence Isoelectric_std', 'Sequence Hyrophobicity_mean','Sequence Hyrophobicity_std', 
            'Sequence Acidity_mean','Sequence Acidity_std','Sequence Aromatic Percentage',
            'Sequence Aromatic Percentage Start', 'Sequence Aromatic Percentage End']

        features = features + ['Amino Acid Percentage_' + str(x) for x in range(22)] #Using 22 since X ignored
        features = features + ['Amino Acid Percentage Start_' + str(x) for x in range(22)]
        features = features + ['Amino Acid Percentage End_' + str(x) for x in range(22)]

        conv_features = []
        for col in self.data.columns:
            if 'k3_' in col or 'k5_' in col or 'k7_' in col:
                conv_features.append(col)
        features = features + conv_features
    
        self.X = self.data[features]
        self.Y = self.data['Label']
        LE = preprocessing.LabelEncoder()
        self.Y = LE.fit_transform(self.Y)
        
    def create_train_test_set(self, ordering, split=0.15):
        data_size = self.Y.shape[0]
        val_indices = ordering[:round(split*data_size)]
        train_indices = ordering[round(split*data_size):]
        train_x = self.X.iloc[train_indices]
        train_y = self.Y[train_indices]
        val_x = self.X.iloc[val_indices]
        val_y = self.Y[val_indices]
        return train_x, train_y, val_x, val_y
