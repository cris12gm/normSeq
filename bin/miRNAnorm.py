import sys,os
import numpy as np
import pandas as pd

def readTxt(fileInput,fileOutput):
    cabecera = open(fileInput).readline().split("\t")[0]
    df = pd.read_table(sys.argv[1])
    df = df.set_index(cabecera) ##Esto hay que generalizarlo para que coja lo que haya en la primera columna

    #Normalize RPM
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    normalized = np.divide(df, sums)*1000000

    normalized.to_csv(fileOutput,sep="\t")

readTxt(sys.argv[1],sys.argv[2])
