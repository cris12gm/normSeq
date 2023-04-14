import os,sys
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler


from config import METHODS,R_PATH,R_SCRIPTS_PATH
import subprocess

def createplots(infile,df,method,jobDir,annotation,annotation_df,cmds):
    outDir = os.path.join(jobDir,"graphs")

    if not os.path.exists(outDir):
        os.mkdir(outDir)

    # PCA
    outfile = os.path.join(outDir,"pca_"+method+".html")
    outfileImage = os.path.join(outDir,"pca_"+method+".png")
    pca(df,annotation_df,outfile,outfileImage)

    # Heatmap
    outfile = os.path.join(outDir,"heatmap_"+method+".png")
    outfileImage = os.path.join(outDir,"heatmap_"+method+".png")
    cmd_heatmap = heatmap(infile,outfile,annotation)
    cmds.append(cmd_heatmap)

    return(cmds)
    
def heatmap(infile,outfile,annotation):

    #Execute in R
    cmd_heatmap = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"heatmap.R "+infile+" "+annotation+" "+outfile
    return(cmd_heatmap)

def pca(df,annotation_df,outfile,outfileImage):

    dfT = df.T
    norm_X = MinMaxScaler().fit_transform(dfT)
    pca = PCA(n_components=2)
    components = pca.fit_transform(norm_X)
    variance = pca.explained_variance_ratio_*100
    total_var = pca.explained_variance_ratio_.sum() * 100


    fig = px.scatter(components, x=0, y=1, hover_name=list(annotation_df.index), color=list(annotation_df.group),
        title=f'Total Explained Variance: {total_var:.2f}%')

    fig.update_layout(
    xaxis_title="PC1 ("+str(round(variance[0],2))+"%)",
    yaxis_title="PC2 ("+str(round(variance[1],2))+"%)",
    legend_title="Group",
    font=dict(
        size=18
    )
    )

    fig.update_traces(marker_size=10)
    fig.write_image(outfileImage)

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()


    df = px.data.iris()

    pca = PCA(n_components=3)
    components = pca.fit_transform(norm_X)

    total_var = pca.explained_variance_ratio_.sum() * 100

    fig3D = px.scatter_3d(
        components, x=0, y=1, z=2, color=list(annotation_df.group),
        title=f'Total Explained Variance: {total_var:.2f}%',
        labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
    )
    fig3D.update_layout(
    legend_title="Group"
    )
    
    plotCode = plot(fig3D, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile.replace(".html","_3D.html"),'a')
    outfile_W.write(plotCode)
    outfile_W.close()



def pca_batch(df,annotation_df,outfile,outfileImage):

    dfT = df.T
    norm_X = MinMaxScaler().fit_transform(dfT)
    pca = PCA(n_components=2)
    components = pca.fit_transform(norm_X)
    variance = pca.explained_variance_ratio_*100
    total_var = pca.explained_variance_ratio_.sum() * 100


    fig = px.scatter(components, x=0, y=1, hover_name=list(annotation_df.index), color=list(annotation_df.batchEffect),
        title=f'Total Explained Variance: {total_var:.2f}%')

    fig.update_layout(
    xaxis_title="PC1 ("+str(round(variance[0],2))+"%)",
    yaxis_title="PC2 ("+str(round(variance[1],2))+"%)",
    legend_title="Group",
    font=dict(
        size=18
    )
    )

    fig.update_traces(marker_size=10)
    fig.write_image(outfileImage)

    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile,'a')
    outfile_W.write(plotCode)
    outfile_W.close()


    df = px.data.iris()

    pca = PCA(n_components=3)
    components = pca.fit_transform(norm_X)

    total_var = pca.explained_variance_ratio_.sum() * 100

    fig3D = px.scatter_3d(
        components, x=0, y=1, z=2, color=list(annotation_df.group),
        title=f'Total Explained Variance: {total_var:.2f}%',
        labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
    )
    fig3D.update_layout(
    legend_title="Group"
    )
    
    plotCode = plot(fig3D, show_link=False, auto_open=False, output_type = 'div')
    outfile_W = open(outfile.replace(".html","_3D.html"),'a')
    outfile_W.write(plotCode)
    outfile_W.close()

