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

def createsummary(infile,df,method,jobDir,annotation):
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
    distribution(df,outfile,outfileImage,title)

def distribution(df,outfile,outfileImage,title):

    fig = make_subplots(rows=1, cols=1)
    hist_data = []
    group_labels = []
    for col in df.columns:
        hist_data_this = np.log2(df[col]+1).values.tolist()
        hist_data.append(hist_data_this)
        group_labels.append(col)
    fig = ff.create_distplot(hist_data,group_labels,show_hist=False,show_rug=False)
    fig.update_layout(title_text=title)
    fig.write_image(outfileImage,scale=3,height=400)
    
    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()
