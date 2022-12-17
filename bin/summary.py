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
import subprocess

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
    distribution(df,outfile,outfileImage)



def distribution(df,outfile,outfileImage):

    fig = make_subplots(rows=1, cols=1)
    for col in df.columns:
        trace0 = go.Histogram(x=np.log2(df[col]+1),legendgroup=col)
        fig.append_trace(trace0, 1, 1)
    fig.show()
    sys.exit(1)
    fig = px.histogram(df, x="Control_M",color="Sample")

    fig.write_image(outfileImage)

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()
