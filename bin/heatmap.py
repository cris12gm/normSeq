import plotly.express as px
from plotly.offline import plot
import plotly.graph_objects as go

import pandas as pd
import numpy as np
import sys

pd.options.plotting.backend = "plotly"

infile = open(sys.argv[1],'r')

df = pd.read_table(infile,sep="\t",index_col="name")
rpm = df.values.T
sumrpm =  np.array(pd.DataFrame(np.sum(rpm,axis=0)).T)

normalized = np.divide(rpm, sumrpm).T
original = rpm.T

samples = df.columns.tolist()
mirnas = df.index.tolist()


fig = go.Figure(data=go.Heatmap(
    z = original,
    x = samples,
    y = mirnas,
    text = original,
    hoverinfo = 'x+y+text',
    colorscale='Viridis')
)


#fig = px.imshow(df, color_continuous_scale='RdBu_r', origin='lower')


plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
outfile =  open(sys.argv[2],'a')
outfile.write(plotCode)
outfile.close()