import os,sys
import pandas as pd
import numpy as np
import subprocess
from plots import pca
from config import R_PATH
from config import R_SCRIPTS_PATH
import itertools

def processAnnotation(infile):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df.rename(columns = {cabecera:'sample'}, inplace = True)
    df = df.set_index(cabecera)
    df = df.dropna()
    return(df)

def edgeR(infile,method,annotation,FDR,jobDir):
    #Get groups from annotation

    if not os.path.exists(os.path.join(jobDir,"DE")):
        os.mkdir(os.path.join(jobDir,"DE"))

    annotationFile = processAnnotation(annotation)
    groups = annotationFile["group"].values.tolist()
    diffGroups = list(set(groups))

    for subset in itertools.combinations(diffGroups, 2):
        group1 = subset[0]
        group2 = subset[1]
        output = os.path.join(jobDir,"DE","edgeR_"+group1+"_"+group2+".txt")
        subprocess.call (R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_de.R "+infile+" "+method+" "+annotation+" "+FDR+" "+group1+" "+group2+" "+output,shell=True)

