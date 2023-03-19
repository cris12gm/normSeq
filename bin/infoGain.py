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
    infoGain = mutual_info_classif(data.reshape(-1, 1), target,discrete_features=False,random_state=1,n_neighbors=3)[0]
    if infoGain>1:
        infoGain = 1
    return infoGain

def calculate_infoGain(df,annotation_df,group):

    annotation_group = annotation_df[(annotation_df['group'] == group)]
    samplesGroup = list(annotation_group.index)
    dfGroup = df[samplesGroup]
  
    annotation_other = annotation_df[(annotation_df['group'] != group)]
    annotation_other = annotation_other.drop('group',axis=1)
    annotation_other['group']="Other"
    samplesOther = list(annotation_other.index)
    dfOther = df[samplesOther]
    annotationC = pd.concat([annotation_group,annotation_other])
    dfC = pd.merge(dfGroup,dfOther,left_index=True,right_index=True)
    merged = (pd.merge(dfC.T,annotationC,left_index=True,right_index=True)).T
    labels = np.asarray((merged.loc['group']).values)
    merged = np.asarray(merged.drop("group").T)

    # Calculate information gain for each feature in parallel using joblib
    results = Parallel(n_jobs=12)(delayed(calc_information_gain)(merged[:, i], labels) for i in range(merged.shape[1]))

    return(results)

def plotInfo(df,outfileImage,outfile,titleThis,title):

    
    fig = px.box(df,title=titleThis,notched=True)
    fig.update_layout(
    xaxis_title="Method",
    yaxis_title=title,yaxis_range=[0,1])
    
    fig.write_image(outfileImage)

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()
