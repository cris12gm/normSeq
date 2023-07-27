import numpy as np
from sklearn.datasets import load_iris
from sklearn.feature_selection import mutual_info_classif
from joblib import Parallel, delayed
import sys
import pandas as pd
import plotly.express as px
from plotly.offline import plot
import multiprocessing as mp
import plotly.graph_objects as go

def calc_information_gain(data, target):
    # Calculate information gain for a single feature using mutual_info_classif
    errorInfo = False
    try:
        infoGain = mutual_info_classif(data.reshape(-1, 1), target,discrete_features=False,random_state=1,n_neighbors=3)[0]
    except:
        errorInfo = True
    if infoGain>1:
        infoGain = 1
    return infoGain,errorInfo

def calculate_infoGain_group(df,annotation_df,group):
    try:
        annotation_df = annotation_df.drop("replicate",axis=1)
    except:
        pass

    annotation_group = annotation_df[(annotation_df['group'] == group)]
    samplesGroup = list(annotation_group.index)
    dfGroup = df[samplesGroup]

    annotation_other = annotation_df[(annotation_df['group'] != group)]
    annotation_other = annotation_other.drop('group',axis=1)
    annotation_other['group']="Other"
    samplesOther = list(annotation_other.index)
    if (len(samplesGroup) + len(samplesOther)) <3:
        results=False
    else:
        dfOther = df[samplesOther]
        annotationC = pd.concat([annotation_group,annotation_other])
        dfC = pd.merge(dfGroup,dfOther,left_index=True,right_index=True)
        merged = (pd.merge(dfC.T,annotationC,left_index=True,right_index=True)).T
        labels = ((merged.loc['group']).values)
        merged = merged.drop("group",axis=0).T
        merged = (merged).to_numpy()
        # Calculate information gain for each feature in parallel using joblib
        results,errorInfo = Parallel(n_jobs=2)(delayed(calc_information_gain)(merged[:, i], labels) for i in range(merged.shape[1]))
        if errorInfo:
            results="Error"
    return(results)


def calculate_infoGain_pairwise(df,annotation_df,combination):

    try:
        annotation_df = annotation_df.drop("replicate",axis=1)
    except:
        pass
    group1 = combination[0]
    group2 = combination[1]

    filtered_Annotation = annotation_df[(annotation_df['group'] == group1) | (annotation_df['group'] == group2)]
    selectedSamples = list(filtered_Annotation.index)
    if len(selectedSamples)>2:
        df2 = df[selectedSamples]
        merged = (pd.merge(df2.T,filtered_Annotation,left_index=True,right_index=True)).T
        labels = np.asarray((merged.loc['group']).values)
        merged = np.asarray(merged.drop("group").T)
        # Calculate information gain for each feature in parallel using joblib
        results,errorInfo = Parallel(n_jobs=12)(delayed(calc_information_gain)(merged[:, i], labels) for i in range(merged.shape[1]))
        if errorInfo:
            results="Error"
    else:
        results = False
    return(results)


def plotInfo(df,outfileImage,outfile,titleThis,title):

    
    fig = px.box(df,title=titleThis,notched=True)
    fig.update_layout(
    xaxis_title="Method",
    yaxis_title=title,
    font=dict(
        size=16
    ))
    
    fig.write_image(outfileImage)

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()

# def top10Info(infoDf,normalized,methods,annotation_df,outfileImage,outfile):
#     thisDf = infoDf

#     thisDf["max"] = thisDf.max(axis=1)
#     thisDf["min"] = thisDf.min(axis=1)
#     thisDf["diff"] = thisDf["max"] - thisDf["min"]
#     top10 = list((thisDf.sort_values('diff').head(10)).index)

#     exprPlot = []
#     for method in methods:
#         normDf = normalized[method][0]
#         result = normDf.loc[top10].T
#         result = (pd.merge(result,annotation_df,left_index=True,right_index=True))
#         result = result.set_index('group')
#         result = pd.melt(result)
#         result["method"] = method
#         exprPlot.append(result)
#     expression = pd.concat(exprPlot)

#     print(expression)
#     fig = go.Figure()
#     fig = px.box(expression,color='variable')
#     # labels={
#     #                  "group": "Group",
#     #                  "value":"Expression",
#     #                  "variable":""
#     #              },)
    
#     fig.show()


#     # sys.exit(1)
#     # fig.write_image(outfileImage,scale=3,height=400)
    
#     # plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
#     # outfile_W = open(outfile,'a')
#     # outfile_W.write(plotCode)
#     # outfile_W.close()