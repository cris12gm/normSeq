import numpy as np
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go

def heatmap(infile,outfile):
    pd.options.plotting.backend = "plotly"
    df = pd.read_table(infile,sep="\t",index_col="name")
    rpm = df.values.T
    sumrpm =  np.array(pd.DataFrame(np.sum(rpm,axis=0)).T)

    normalized = np.divide(rpm, sumrpm).T
    original = rpm.T

    samples = df.columns.tolist()
    mirnas = df.index.tolist()


    fig = go.Figure(data=go.Heatmap(
        z = normalized,
        x = samples,
        y = mirnas,
        text = original,
        hoverinfo = 'x+y+text',
        colorscale='Viridis')
    )

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()