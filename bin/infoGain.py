import numpy as np
from sklearn.datasets import load_iris
from sklearn.feature_selection import mutual_info_classif
from joblib import Parallel, delayed
import sys
import pandas as pd
import plotly.express as px
from plotly.offline import plot
import multiprocessing as mp

def calc_information_gain(data, target):
    # Calculate information gain for a single feature using mutual_info_classif
    infoGain = mutual_info_classif(data.reshape(-1, 1), target)[0]
    return infoGain

def calculate_infoGain(df,annotation_df,combination):
    group1 = combination[0]
    group2 = combination[1]

    filtered_Annotation = annotation_df[(annotation_df['group'] == group1) | (annotation_df['group'] == group2)]
    selectedSamples = list(filtered_Annotation.index)
    df2 = df[selectedSamples]

    merged = (pd.merge(df2.T,filtered_Annotation,left_index=True,right_index=True)).T
    labels = np.asarray((merged.loc['group']).values)
    merged = np.asarray(merged.drop("group").T)

    # Calculate information gain for each feature in parallel using joblib
    results = Parallel(n_jobs=12)(delayed(calc_information_gain)(merged[:, i], labels) for i in range(merged.shape[1]))

    return(results)

def plotInfo(df,outfileImage,outfile,titleThis):

    fig = px.box(df,title=titleThis)
    fig.update_layout(
    xaxis_title="Method",
    yaxis_title="Information Gain per miRNA")
    
    fig.write_image(outfileImage)

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()
