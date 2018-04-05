import numpy as np
import pandas as pd
from Bio import SeqIO
from random import shuffle


class DataLoader(object):


    def __init__(self, path='./'):
        self.path = path
        self.labels = []
        self.sequences = []
        self.data = None
        
    def load_files(self, files, labels):
        file_index = 0
        for file in files:
            self.load_file(file, labels[file_index])
            file_index += 1
        
    def load_file(self, file, label):
        print("Loading file: {}...".format(self.path+file))
        with open(self.path+file, 'r' ) as f:
            content = SeqIO.parse(f, 'fasta')
            for record in content:
                seq = str(record.seq)
                self.labels.append(label)
                self.sequences.append(seq)
        print("... Success")
    
    def convert_to_df(self, verbose=True):
        print("Converting data to pandas DataFrame")
        self.data = pd.DataFrame({'Label': self.labels, 'Sequence': self.sequences})
        print("Success")
        if verbose:
            print("Example data:")
            print(self.data.head())
    
        
        
    
        
    