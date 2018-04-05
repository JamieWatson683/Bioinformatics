import pandas as pd
import numpy as np

#Helper Functions

def data_viewer(data):
    shape = data.shape
    max_len = data.Sequence.map(len).max()
    avg_len = data.Sequence.map(len).mean()
    std_len = data.Sequence.map(len).std()
    min_len = data.Sequence.map(len).min()
    label_counts = data.Label.value_counts()
    
    print("Shape: {}".format(shape))
    print("Max Seq Length: {}".format(max_len))
    print("Min Seq Length: {}".format(min_len))
    print("Avg Seq Length: {}".format(avg_len))
    print("Standard Deviation: {}".format(std_len))
    print("Sequence Class Counts -->")
    print(label_counts)
    
def split_train_validation(df, split=0.2):
    #Randomly Shuffle
    df = df.sample(frac=1.0)
    #Split into training and validation set
    rows = df.shape[0]
    validation = df[:round(rows*split)]
    train = df[round(rows*split):]
    return train, validation