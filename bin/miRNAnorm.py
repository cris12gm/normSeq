import sys,os
import numpy as np
import pandas as pd

def readTxt(infile,outfile):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df = df.set_index(cabecera) ##Esto hay que generalizarlo para que coja lo que haya en la primera columna

    #Normalize RPM
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    normalized = np.divide(df, sums)*1000000

    normalized.to_csv(outfile,sep="\t")
