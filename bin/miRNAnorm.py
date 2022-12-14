import os,sys
import pandas as pd
import numpy as np
import subprocess

from config import R_PATH

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

    return normalized

def tc(df,outfile):
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    mean = np.mean(sums)

    normalized = np.divide(df,mean)
    normalized.to_csv(outfile,sep="\t")

    return normalized

def uq(infile,outfile):
    #Execute in R
    subprocess.call (R_PATH+" --vanilla bin/R/edgeR_normalization.R "+infile+" UQ "+outfile, shell=True)
    outdf = processInput(outfile)
    print (outdf)
    return outdf

def tmm(infile,outfile):
    #Execute in R
    subprocess.call (R_PATH+" --vanilla bin/R/edgeR_normalization.R "+infile+" TMM "+outfile, shell=True)
    outdf = processInput(outfile)
    return outdf

def rle(infile,outfile):
    #Execute in R
    subprocess.call (R_PATH+" --vanilla bin/R/edgeR_normalization.R "+infile+" RLE "+outfile, shell=True)
    outdf = processInput(outfile)
    return outdf

def med(df,outfile):
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    median = np.median(sums)

    normalized = np.divide(df,median)
    normalized.to_csv(outfile,sep="\t")

    return normalized

def deseq(infile,outfile):
    #Execute in R
    subprocess.call (R_PATH+" --vanilla bin/R/deseq_normalization.R "+infile+" "+outfile, shell=True)
    outdf = processInput(outfile)
    return outdf


def qn(infile,outfile):
    #Execute in R
    subprocess.call (R_PATH+" --vanilla bin/R/quantile_normalization.R "+infile+" "+outfile, shell=True)
    outdf = processInput(outfile)
    return outdf

def ruv(infile,outfile):
    #Execute in R
    subprocess.call (R_PATH+" --vanilla bin/R/ruv_normalization.R "+infile+" "+outfile, shell=True)
    outdf = processInput(outfile)
    return outdf

def norm(infile,df,method,jobDir):
    outDir = os.path.join(jobDir,"normalized")
    print (method)
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    if method == "CPM":
        outfile = os.path.join(outDir,"matrix_CPM.txt")
        outdf = cpm(df,outfile)
    elif method =="TC":
        outfile = os.path.join(outDir,"matrix_TC.txt")
        outdf = tc(df,outfile)
    elif method == "UQ":
        outfile = os.path.join(outDir,"matrix_UQ.txt")
        print (outfile)
        outdf = uq(infile,outfile)
        
    elif method == "TMM":
        outfile = os.path.join(outDir,"matrix_TMM.txt")
        outdf = tmm(infile,outfile)
    elif method == "RLE":
        outfile = os.path.join(outDir,"matrix_RLE.txt")
        outdf = rle(infile,outfile)
    elif method == 'NN':
        outfile = os.path.join(jobDir,"matrix.txt")
        outdf = df
    elif method == 'Med':
        outfile = os.path.join(outDir,"matrix_Med.txt")
        outdf = med(df,outfile)
    elif method == "DESEQ":
        outfile  = os.path.join(outDir,"matrix_DESEQ.txt")
        outdf = deseq(infile,outfile)
    elif method == "QN":
        outfile = os.path.join(outDir,"matrix_QN.txt")
        outdf = qn(infile,outfile)

    return outdf,outfile
