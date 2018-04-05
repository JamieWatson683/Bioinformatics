import pandas as pd
import numpy as np


LOOKUP = {'A': 0, 'R': 1, 'N':2, 'D': 3, 'C': 4, 'E':5, 'Q': 6, 'G': 7, 'H': 8, 'I':9, 'L': 10, 'K':11, 'M': 12,
            'F':13, 'P':14, 'S':15, 'T': 16, 'W':17, 'Y':18, 'V':19, 'U': 20,'B':21}

#AMINO_VECTOR = np.array([[71, 6, .806, 2.34], # [Weight, Isoelectrc point, Hydrophobicity, Acidity]
#                          [156, 10.8, 0,2.17],
#                          [114,5.4,.448,2.02],
#                          [115,3,.417,1.88],
#                          [103,5,.721,1.96],
#                          [129,3.2,.458,2.19],
#                          [128,5.7,.43,2.17],
#                          [57,6,.77,2.34],
#                          [137,7.6,.548,1.82],
#                          [113,6,1,2.36],
#                         [113,6,.918,2.36],
#                          [128,9.7,.263,2.18],
#                          [131,5.7,.811,2.28],
#                          [147,5.5,.951,1.83],
#                          [97,6.3,.678,1.99],
#                          [87,5.7,.601,2.21],
#                          [101,5.6,.634,2.09],
#                          [186,5.9,.854,2.83],
#                          [163,5.7,.714,2.20],
#                          [99,6,.923,2.32],
#                          [121,5.7,.64725,2.177], #U-> defined weight and charge, mean hydro and acidity
#                          [119,6,.64725,2.177], #X-> mean of all values
#                          [114.5,4.2,0.4352,1.95]]) #B-> mean of N,D


# [Weight, Isoelectrc point, Hydrophobicity, pKa, aromatic (0,1)]
AMINO_VECTOR = np.array([[71, 6, 1.8, 2.34, 0], 
                          [156, 10.8, -4.5, 2.17, 0],
                          [114,5.4,-3.5,2.02, 0],
                          [115,3,-3.5,1.88, 0],
                          [103,5,2.5,1.96, 0],
                          [129,3.2,-3.5,2.19, 0],
                          [128,5.7,-3.5,2.17, 0],
                          [57,6,-.4,2.34, 0],
                          [137,7.6,-3.2,1.82, 1],
                          [113,6,4.5,2.36, 0],
                          [113,6,3.8,2.36, 0],
                          [128,9.7,-3.9,2.18, 0],
                          [131,5.7,1.9,2.28, 0],
                          [147,5.5,2.8,1.83, 1],
                          [97,6.3,-1.6,1.99, 0],
                          [87,5.7,-.8,2.21, 0],
                          [101,5.6,-.7,2.09, 0],
                          [186,5.9,-.9,2.83, 1],
                          [163,5.7,-1.3,2.20, 1],
                          [99,6,4.2,2.32, 0],
                          [121,5.7,2.5,2.177, 0], #U-> defined weight and charge, same hydro as cysteine
                          [114.5,4.2,-3.5,1.95, 0]]) #B-> mean of N,D
    
#Feature Selection Functions

def sequence_length(df):
    print("Adding Sequence Length...")
    seq_lens = df.Sequence.map(len)
    mean_len = seq_lens.mean()
    std_len = seq_lens.std()
    seq_lens_norm = (seq_lens - mean_len) / std_len
    df['Sequence Length'] = seq_lens
    df['Sequence Length Norm'] = seq_lens_norm
    print("Success")
    
def numerate_sequence_and_attributes(df):
    """NOTE: Skips 'X' occurances"""
    print("Numerating Sequences...")
    print("Getting Attributes...")
    num_seq = []
    df_rows = df.shape[0]
    weight = []
    isoelectric = []
    hydrophobicity = []
    acidity = []
    aromatic = []
    for row in range(df_rows):
        seq = []
        weight_ = []
        isoelectric_ = []
        hydrophobicity_ = []
        acidity_ = []
        aromatic_ = []
        for amino in df.Sequence[row]:
            if amino != 'X':
                seq.append(LOOKUP[amino])
                weight_.append(AMINO_VECTOR[seq[-1],0])
                isoelectric_.append(AMINO_VECTOR[seq[-1],1])
                hydrophobicity_.append(AMINO_VECTOR[seq[-1],2])
                acidity_.append(AMINO_VECTOR[seq[-1],3])
                aromatic_.append(AMINO_VECTOR[seq[-1],4])
        num_seq.append(seq)
        weight.append(weight_)
        isoelectric.append(isoelectric_)
        hydrophobicity.append(hydrophobicity_)
        acidity.append(acidity_)
        aromatic.append(aromatic_)
    df['Sequence Numerated'] = num_seq
    df['Sequence Weights'] = weight
    df['Sequence Isoelectric'] = isoelectric
    df['Sequence Hyrophobicity'] = hydrophobicity
    df['Sequence Acidity'] = acidity
    df['Sequence Aromatic'] = aromatic
    print("Success")
            
def amino_acid_count(df, first_size=50, last_size=50):
    print("Adding Amino Acid Counts...")
    df_rows = df.shape[0]
    amino_acids = len(LOOKUP.keys())
    counts = []
    first_counts = []
    last_counts = []
    for row in range(df_rows):
        seq_counts = np.zeros(amino_acids)
        seq_counts_first = np.zeros(amino_acids)
        seq_counts_last = np.zeros(amino_acids)
        for index in range(amino_acids):
            seq_counts[index] = df['Sequence Numerated'][row].count(index)
            seq_counts_first[index] = df['Sequence Numerated'][row][:first_size].count(index)
            seq_counts_last[index] = df['Sequence Numerated'][row][-last_size:].count(index)
        counts.append(seq_counts)
        first_counts.append(seq_counts_first)
        last_counts.append(seq_counts_last)
    df['Amino Acid Count'] = counts
    df['Amino Acid Count Start'] = first_counts
    df['Amino Acid Count End'] = last_counts
    print("Success")
    
    
def amino_acid_percentage(df, first_size=50, last_size=50):
    print("Adding Percentages...")
    percentage = df['Amino Acid Count'] / df['Sequence Length']
    percentage_first = df['Amino Acid Count Start'] / 50
    percentage_last = df['Amino Acid Count End'] / 50
    df['Amino Acid Percentage'] = percentage
    df['Amino Acid Percentage Start'] = percentage_first
    df['Amino Acid Percentage End'] = percentage_last
    print("Success")    

def explode_feature(df, feature):
    print("Exploding feature: {}".format(feature))
    num_rows = df.shape[0]
    aminos = len(LOOKUP.keys())
    array = np.zeros((num_rows, aminos))
    for row in range(num_rows):
        array[row] = df[feature].iloc[row]
    array = pd.DataFrame(array)
    array.columns = [feature + '_' + str(x) for x in range(aminos)]
    print("Success")
    return pd.concat([df, array], axis=1)

def expand_attributes(df):
    expand_feature_list(df, 'Sequence Weights')
    expand_feature_list(df, 'Sequence Isoelectric')
    expand_feature_list(df, 'Sequence Hyrophobicity')
    expand_feature_list(df, 'Sequence Acidity')

def expand_feature_list(df, feature):
    print("Expanding Attribute: {}...".format(feature))
    df[feature+"_mean"] = df[feature].map(np.mean)
    df[feature+"_std"] = df[feature].map(np.std)
    print("Convolving...")
    df[feature+"_k3"] = df[feature].map(convolve_feature_k3)
    df[feature+"_k3_mean"] = df[feature+"_k3"].map(np.mean)
    df[feature+"_k3_std"] = df[feature+"_k3"].map(np.std)
    df[feature+"_k3_max"] = df[feature+"_k3"].map(np.max)
    df[feature+"_k3_min"] = df[feature+"_k3"].map(np.min)
    df[feature+"_k3_far_abv"] = df[feature+"_k3"].map(far_from_mean_above)
    df[feature+"_k3_far_blw"] = df[feature+"_k3"].map(far_from_mean_below)
    print("k3 Done")
    df[feature+"_k5"] = df[feature].map(convolve_feature_k5)
    df[feature+"_k5_mean"] = df[feature+"_k5"].map(np.mean)
    df[feature+"_k5_std"] = df[feature+"_k5"].map(np.std)
    df[feature+"_k5_max"] = df[feature+"_k5"].map(np.max)
    df[feature+"_k5_min"] = df[feature+"_k5"].map(np.min)
    df[feature+"_k5_far_abv"] = df[feature+"_k5"].map(far_from_mean_above)
    df[feature+"_k5_far_blw"] = df[feature+"_k5"].map(far_from_mean_below)
    print("k5 Done")
    df[feature+"_k7"] = df[feature].map(convolve_feature_k7)
    df[feature+"_k7_mean"] = df[feature+"_k7"].map(np.mean)
    df[feature+"_k7_std"] = df[feature+"_k7"].map(np.std)
    df[feature+"_k7_max"] = df[feature+"_k7"].map(np.max)
    df[feature+"_k7_min"] = df[feature+"_k7"].map(np.min)
    df[feature+"_k7_far_abv"] = df[feature+"_k7"].map(far_from_mean_above)
    df[feature+"_k7_far_blw"] = df[feature+"_k7"].map(far_from_mean_below)
    print("k7 Done")
    
    
def convolve_feature_k3(feature_data):
    kernel = np.ones(3)
    feature_data = np.array(feature_data)
    data_len = feature_data.shape[0]
    result = np.zeros(data_len-kernel.shape[0]+1)
    for i in range(0,data_len-kernel.shape[0]+1):
        result[i] = np.dot(feature_data[i:i+kernel.shape[0]], kernel)
    return result

def convolve_feature_k5(feature_data):
    kernel = np.ones(5)
    feature_data = np.array(feature_data)
    data_len = feature_data.shape[0]
    result = np.zeros(data_len-kernel.shape[0]+1)
    for i in range(0,data_len-kernel.shape[0]+1):
        result[i] = np.dot(feature_data[i:i+kernel.shape[0]], kernel)
    return result

def convolve_feature_k7(feature_data):
    kernel = np.ones(7)
    feature_data = np.array(feature_data)
    data_len = feature_data.shape[0]
    result = np.zeros(data_len-kernel.shape[0]+1)
    for i in range(0,data_len-kernel.shape[0]+1):
        result[i] = np.dot(feature_data[i:i+kernel.shape[0]], kernel)
    return result

def far_from_mean_above(feature):
    size = len(feature)
    mean = np.mean(feature)
    std = np.std(feature)
    count = np.sum(feature > (mean + 2*std))
    return count / size

def far_from_mean_below(feature):
    size = len(feature)
    mean = np.mean(feature)
    std = np.std(feature)
    count = np.sum(feature < (mean - 2*std))
    return count / size

def aromatic_percentage(df):
    print("Extracting aromatic features...")
    global_percentage = df['Sequence Aromatic'].map(np.sum) / df['Sequence Length']
    start_percentage = df['Sequence Aromatic'].map(aromatic_start)
    end_percentage = df['Sequence Aromatic'].map(aromatic_end)
    df['Sequence Aromatic Percentage'] = global_percentage
    df['Sequence Aromatic Percentage Start'] = start_percentage
    df['Sequence Aromatic Percentage End'] = end_percentage
    print("Done")
    
def aromatic_start(aromatic_data):
    return np.sum(aromatic_data[0:50]) / 50

def aromatic_end(aromatic_data):
    return np.sum(aromatic_data[-50:]) / 50