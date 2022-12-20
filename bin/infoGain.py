from scipy.stats import entropy
import pandas as pd
import numpy as np
import sys

def compute_impurity(feature, impurity_criterion):
    """
    This function calculates impurity of a feature.
    Supported impurity criteria: 'entropy', 'gini'
    input: feature (this needs to be a Pandas series)
    output: feature impurity
    """
    probs = feature.value_counts(normalize=True)
    
    if impurity_criterion == 'entropy':
        impurity = -1 * np.sum(np.log2(probs) * probs)
    elif impurity_criterion == 'gini':
        impurity = 1 - np.sum(np.square(probs))
    else:
        raise ValueError('Unknown impurity criterion')
        
    return(round(impurity, 3))

def comp_feature_information_gain(df, target, descriptive_feature, split_criterion):
    """
    This function calculates information gain for splitting on 
    a particular descriptive feature for a given dataset
    and a given impurity criteria.
    Supported split criterion: 'entropy', 'gini'
    """
    
    # print('target feature:', target)
    # print('descriptive_feature:', descriptive_feature)
    # print('split criterion:', split_criterion)
            
    target_entropy = compute_impurity(df[target], split_criterion)

    # we define two lists below:
    # entropy_list to store the entropy of each partition
    # weight_list to store the relative number of observations in each partition
    entropy_list = list()
    weight_list = list()
    
    # loop over each level of the descriptive feature
    # to partition the dataset with respect to that level
    # and compute the entropy and the weight of the level's partition
    for level in df[descriptive_feature].unique():
        df_feature_level = df[df[descriptive_feature] == level]
        entropy_level = compute_impurity(df_feature_level[target], split_criterion)
        entropy_list.append(round(entropy_level, 3))
        weight_level = len(df_feature_level) / len(df)
        weight_list.append(round(weight_level, 3))

    # print('impurity of partitions:', entropy_list)
    # print('weights of partitions:', weight_list)

    feature_remaining_impurity = np.sum(np.array(entropy_list) * np.array(weight_list))
    # print('remaining impurity:', feature_remaining_impurity)
    
    information_gain = target_entropy - feature_remaining_impurity
    # print('information gain:', information_gain)
    
    # print('====================')

    return(information_gain)


infile = 'C:/Users/Cris/Dropbox/TRABAJO/NL/normSeq/Website/normSeq/uploads/HYLM9GQEDKRG69K/matrix.txt'
cabecera = open(infile).readline().split("\t")[0]
df = pd.read_table(infile)
df.rename(columns = {cabecera:'name'}, inplace = True)
df = df.set_index(cabecera)
df = df.dropna()
df = df.T


#Read annotation
annotationFile = "C:/Users/Cris/Dropbox/TRABAJO/NL/normSeq/Website/normSeq/uploads/HYLM9GQEDKRG69K/annotation.txt"
cabecera = open(annotationFile).readline().split("\t")[0]
ann = pd.read_table(annotationFile)
ann.rename(columns = {cabecera:'sample'}, inplace = True)
ann = ann.set_index(cabecera)


together = pd.concat([df, ann], axis=1, join="inner")


split_criterion = 'entropy'
infoGain = {}
for feature in together.drop(columns='group').columns:
    feature_info_gain = comp_feature_information_gain(together, 'group', feature, split_criterion)
    infoGain[feature]=feature_info_gain

data = pd.DataFrame.from_dict(infoGain,orient='index')
joined = pd.concat([data, together.T.drop("group")], axis=1, join="inner")
print (joined)