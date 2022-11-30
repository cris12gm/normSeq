import os
import pandas as pd
import numpy as np

def processInput(infile):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df = df.set_index(cabecera)
    df = df.dropna()

    return(df)

def cpm(df,outfile):

    #Normalize CPM
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    normalized = np.divide(df, sums)*1000000

    normalized.to_csv(outfile,sep="\t")

def tc(df,outfile):
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    mean = np.mean(sums)

    normalized = np.divide(df,mean)
    normalized.to_csv(outfile,sep="\t")

def norm(df,method,jobDir):
    outDir = os.path.join(jobDir,"normalized")

    if not os.path.exists(outDir):
        os.mkdir(outDir)

    if method == "CPM":
        outfile = os.path.join(outDir,"matrix_CPM.txt")
        cpm(df,outfile)
    elif method =="TC":
        outfile = os.path.join(outDir,"matrix_TC.txt")
        tc(df,outfile)