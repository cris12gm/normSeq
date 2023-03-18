import pandas as pd
import numpy as np
import plotly.express as px
from plotly.offline import plot
from sklearn.feature_selection import mutual_info_classif
import sys
import timeit

def calculate_infoGain(df,annotation_df,combination):
    group1 = combination[0]
    group2 = combination[1]

    filtered_Annotation = annotation_df[(annotation_df['group'] == group1) | (annotation_df['group'] == group2)]
    selectedSamples = list(filtered_Annotation.index)
    df2 = df[selectedSamples]
    start = timeit.timeit()
    a = mutual_info_classif(df2.T, np.ravel(filtered_Annotation))
    end = timeit.timeit()
    print(end - start)
    info = df2
    info['infoGain'] = a
    infoThis = info['infoGain']
    return infoThis

def plotInfo(df,outfileImage,outfile,titleThis):
    x = []
    y = []
    for element in df:
        x.append(element)
        y.append(df[element])

    fig = px.box(df,title=titleThis)
    fig.update_layout(
    xaxis_title="Method",
    yaxis_title="Information Gain per miRNA")
    
    fig.write_image(outfileImage)

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()