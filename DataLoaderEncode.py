import numpy as np
import pandas as pd
from Bio import SeqIO


class DataLoader(object):


    def __init__(self, path='./', verbose=True, ignoreX=True):
        self.path = path
        self.verbose = verbose
        self.ignoreX = ignoreX
        self.sequences = []
        self.labels = []
        self.seq_lengths = []
        self.max_length = 0
        self.lookup = {'A': 1, 'R': 2, 'N':3, 'D': 4, 'C': 5, 'E':6, 'Q': 7, 'G': 8, 'H': 9, 'I':10, 'L': 11, 'K':12, 'M': 13,
            'F':14, 'P':15, 'S':16, 'T': 17, 'W':18, 'Y':19, 'V':20, 'U': 21, 'B':22, 'X':23}

    def load_files(self, files, labels):
        print("Loading files...")
        index=0
        for file in files:
            if self.verbose:
                print("File {}: {} with label {}".format(index, file, labels[index]))
            self.parse_file(file, labels[index])
            index += 1
        print("Done")

    def parse_file(self, file, label):
        label = self.onehot_label(label)
        with open(self.path+file, 'r' ) as f:
            content = SeqIO.parse(f, 'fasta')
            for record in content:
                seq = str(record.seq)
                seq, seq_len = self.index_sequence(seq)
                self.labels.append(label)
                self.sequences.append(seq)
                self.seq_lengths.append(seq_len)
                
    def onehot_label(self, label):
        if label == 'cyto':
            return np.array([1,0,0,0])
        elif label == 'mito':
            return np.array([0,1,0,0])
        elif label == 'nucleus':
            return np.array([0,0,1,0])
        elif label == 'secreted':
            return np.array([0,0,0,1])
        
    def index_sequence(self, seq):
        seq_len = len(seq)
        if seq_len > self.max_length:
            self.max_length = seq_len
        seq_indexed = np.zeros(seq_len)
        index = 0
        for letter in seq:
            if letter == 'X' and self.ignoreX:
                seq_len -= 1
            else:
                seq_indexed[index] = self.lookup[letter]
                index += 1
        seq_indexed = seq_indexed[:seq_len] #Cut off trailing zeros from ignoring X (if applicable)
        return seq_indexed, seq_len
    
    def convert_all_to_array(self):
        print("Converting labels, sequences to arrays, and created reversed sequence...")
        data_size = len(self.labels)
        labels = np.zeros((data_size, 4))
        sequences_end_pad = np.zeros((data_size, self.max_length))
        sequences_start_pad = np.zeros((data_size, self.max_length))
        for i in range(0,data_size):
            seq_len = self.seq_lengths[i]
            labels[i,:] = self.labels[i]
            sequences_end_pad[i,:seq_len] = self.sequences[i]
            sequences_start_pad[i,-seq_len:] = self.sequences[i]
        self.labels = labels
        self.sequences = sequences_end_pad
        self.sequences_backwards = np.flip(sequences_start_pad, axis=1)
        self.seq_lengths = np.array(self.seq_lengths)
        print("Done")
        
    def trim_sequences(self, sequence_length):
        print("Trimming forward and backward sequences to length {}".format(sequence_length))
        sequences = self.sequences[:,:sequence_length]
        sequences_backwards = self.sequences_backwards[:,:sequence_length]
        self.seq_lengths[self.seq_lengths > sequence_length] = sequence_length
        self.sequences = sequences
        self.sequences_backwards = sequences_backwards
        print("Done")
        
    def create_train_test_set(self,  ordering, split=0.15):
        """Returns train, val: each a list with labels, sequences, sequence backwards, sequence lengths"""
        print("Splittling data with proportion {} to test set...".format(split))
        rnd_seq = self.sequences[ordering]
        rnd_seq_back = self.sequences_backwards[ordering]
        rnd_labels = self.labels[ordering]
        rnd_lens = self.seq_lengths[ordering]
        #Split into training and validation
        train = []
        val = []
        t, v = self.get_split_data(rnd_labels, split)
        train.append(t)
        val.append(v)
        t, v = self.get_split_data(rnd_seq, split)
        train.append(t)
        val.append(v)
        t, v = self.get_split_data(rnd_seq_back, split)
        train.append(t)
        val.append(v)
        t, v = self.get_split_data(rnd_lens, split)
        train.append(t)
        val.append(v)
        print("Done -> returning train, test")
        return train, val
        
    def get_split_data(self, data, split):
        val = data[:round(split*data.shape[0])]
        train = data[round(split*data.shape[0]):]
        return train, val
    
    def get_random_ordering(self):
        return np.random.permutation(self.labels.shape[0])
        
        
    
        
    