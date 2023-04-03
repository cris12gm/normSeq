import os,sys
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import numpy as np

from plotly.subplots import make_subplots
from config import METHODS,R_PATH,R_SCRIPTS_PATH
import plotly.figure_factory as ff
import seaborn as sns
import matplotlib.pyplot as plt

def createsummary(infile,df,method,jobDir,annotation_df,combinations):
    outDir = os.path.join(jobDir,"graphs","summary")
    plotDir = os.path.join(jobDir,"graphs")
    if not os.path.exists(plotDir):
        os.mkdir(plotDir)
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    # distribution
    outfile = os.path.join(outDir,"distribution_"+method+".html")
    outfileImage = os.path.join(outDir,"distribution_"+method+".png")
    title = METHODS[method]
    distribution(df,annotation_df,outfile,outfileImage,title)

    #Top10
    outfile = os.path.join(outDir,"top10_"+method+".html")
    outfileImage = os.path.join(outDir,"top10_"+method+".png")
    outfileTxt = os.path.join(outDir,"top10_"+method+".txt")
    title = METHODS[method]
    top10(df,outfile,outfileImage,outfileTxt,title,annotation_df)

    #Top10FC
    
    title = METHODS[method]
    top10fc(df,outDir,combinations,method,title,annotation_df)


def distribution(df,annotation_df,outfile,outfileImage,title):

    fig = make_subplots(rows=1, cols=1)
    hist_data = []
    group_labels = []
    colors_labels = []
    dfResults = pd.merge(df.T,annotation_df,left_index=True,right_index=True).T
    groups = dfResults.T['group']
    colors = {}
    named_colorscales = px.colors.qualitative.Plotly
    i = 0
    for element in groups:
        try:
            value = colors[element]
        except:
            colors[element] = named_colorscales[i]
            i = i + 1
    for col in df.columns:
        hist_data_this = np.log10(df[col]+1).values.tolist()
        hist_data.append(hist_data_this)
        name = dfResults[col]['group']
        group_labels.append(name)
        colors_labels.append(colors[name])

    fig = ff.create_distplot(hist_data,group_labels,colors=colors_labels,show_hist=False,show_rug=False)
    names = set()
    fig.for_each_trace(
    lambda trace:
        trace.update(showlegend=False)
        if (trace.name in names) else names.add(trace.name))
    fig.update_layout(title_text=title,xaxis_title="Log10(Expression)",
    yaxis_title="Density",
    font=dict(
        size=16
    ))
    fig.write_image(outfileImage,scale=3,height=400)
    
    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()


def top10(df,outfile,outfileImage,outfileTxt,title,annotation):
    media = df.mean(axis=1)
    top10 = list(media.nlargest(10).index)
    
    # Create a DataFrame object from list
    result = df.loc[top10].T
    result = pd.merge(result, annotation, left_index =True ,right_index=True)

    result.to_csv(outfileTxt,sep="\t")
    
    fig = go.Figure()

    fig = px.box(result,color='group',
    labels={
                     "group": "Group",
                     "value":"Expression",
                     "variable":""
                 },)
    

    fig.write_image(outfileImage,scale=3,height=400)
    
    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()

def top10fc(df,outDir,combinations,method,title,annotation_df):

    for combination in combinations:
        
        titlePlot = title+" "+combination[0]+"-"+combination[1]
        group1 = combination[0]
        group2 = combination[1]
        
        outfile = os.path.join(outDir,"top10FC_"+method+"_"+group1+"-"+group2+".html")
        outfileImage = os.path.join(outDir,"top10FC_"+method+"_"+group1+"-"+group2+".png")
        outFiletxt = os.path.join(outDir,"top10FC_"+method+"_"+group1+"-"+group2+".txt")


        samplesFilterg1 = annotation_df[annotation_df['group'] == group1]
        samplesGroup1 = list(samplesFilterg1.index)
        mean_g1 = df[samplesGroup1].mean(axis=1).to_frame()

        samplesFilterg2 = annotation_df[annotation_df['group'] == group2]
        samplesGroup2 = list(samplesFilterg2.index)
        mean_g2 = df[samplesGroup2].mean(axis=1).to_frame()

        result = pd.merge(mean_g1, mean_g2, left_index =True ,right_index=True)

        result.loc[(result['0_x'] >= 10) | (result['0_y'] >= 10)]

        result['fc'] = np.where(result["0_x"]==0,
                        np.abs(np.log2(result["0_y"])), 
                        np.abs(np.log2(result["0_y"]/result["0_x"])))

        top10 = list(result.nlargest(10,'fc').index)

        # Create a DataFrame object from list
        result_expr = df.loc[top10].T
        result_expr = pd.merge(result_expr, annotation_df, left_index =True ,right_index=True)

        result_expr.to_csv(outFiletxt,sep="\t")
        fig = go.Figure()

        fig = px.box(result_expr,color='group',title=titlePlot,
        labels={
                        "group": "Group",
                        "value":"Expression",
                        "variable":""
                    },)

        fig.write_image(outfileImage,scale=3,height=400)
        
        plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
        outfile_W = open(outfile,'a')
        outfile_W.write(plotCode)
        outfile_W.close()