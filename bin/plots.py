import os
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go

def createplots(df,method,jobDir):
    outDir = os.path.join(jobDir,"graphs")

    if not os.path.exists(outDir):
        os.mkdir(outDir)

    # Heatmap
    outfile = os.path.join(outDir,"heatmap_"+method+".html")
    outfileImage = os.path.join(outDir,"heatmap_"+method+".png")
    heatmap(df,outfile,outfileImage,method)

def heatmap(df,outfile,outfileImage,title):
    pd.options.plotting.backend = "plotly"
    rpm = df.values.T
#    sumrpm =  np.array(pd.DataFrame(np.sum(rpm,axis=0)).T)

#   normalized = np.divide(rpm, sumrpm).T
    original = rpm.T

    samples = df.columns.tolist()
    mirnas = df.index.tolist()


    fig = go.Figure(data=go.Heatmap(
#        z = normalized,
        z = original,
        x = samples,
        y = mirnas,
        text = original,
        hoverinfo = 'x+y+text',
        colorscale='Viridis')
    )

    fig.update_layout(title={'text':title}) 

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()

    fig.write_image(outfileImage)
