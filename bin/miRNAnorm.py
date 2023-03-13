import os,sys
import pandas as pd
import numpy as np
import subprocess

from config import R_PATH
from config import R_SCRIPTS_PATH

def processInput(infile):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df.rename(columns = {cabecera:'name'}, inplace = True)
    df = df.set_index('name')
    df = df.dropna()

    return(df)

def processAnnotation(infile):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df.rename(columns = {cabecera:'sample'}, inplace = True)
    df = df.set_index('sample')
    df.dropna(how='all', axis=1, inplace=True)

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

def uq(infile,outfile,jobDir):
    #Execute in R
    cmd_uq = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_normalization.R "+infile+" UQ "+outfile+" >"+jobDir+"/Log.txt"
    return cmd_uq

def tmm(infile,outfile,jobDir):
    #Execute in R
    cmd_tmm = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_normalization.R "+infile+" TMM "+outfile+" >"+jobDir+"/Log.txt"
    return cmd_tmm

def rle(infile,outfile,jobDir):
    #Execute in R
    cmd_rle = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_normalization.R "+infile+" RLE "+outfile+" >"+jobDir+"/Log.txt"
    return cmd_rle

def med(df,outfile):
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    median = np.median(sums)

    normalized = np.divide(df,median)
    normalized.to_csv(outfile,sep="\t")

    return normalized

def deseq(infile,outfile,jobDir):
    #Execute in R
    cmd_deseq = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"deseq_normalization.R "+infile+" "+outfile+" >"+jobDir+"/Log.txt"
    return cmd_deseq


def qn(infile,outfile,jobDir):
    #Execute in R
    cmd_qn = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"quantile_normalization.R "+infile+" "+outfile+" >"+jobDir+"/Log.txt"
    return cmd_qn

def ruv(infile,outfile,jobDir):
    cmd_ruv = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"ruv_normalization.R "+infile+" "+outfile+" >"+jobDir+"/Log.txt"
    return cmd_ruv

def norm(infile,df,method,jobDir):
    outDir = os.path.join(jobDir,"normalized")
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    if method == "CPM":
        outfile = os.path.join(outDir,"matrix_CPM.txt")
        outdf = cpm(df,outfile)
    elif method =="TC":
        outfile = os.path.join(outDir,"matrix_TC.txt")
        outdf = tc(df,outfile)
    elif method == 'NN':
        outfile = os.path.join(outDir,"matrix_NN.txt") 
        outdf = df
    elif method == 'Med':
        outfile = os.path.join(outDir,"matrix_Med.txt")
        outdf = med(df,outfile)

    return outdf,outfile


def norm_r(infile,method,jobDir,cmds):
    outDir = os.path.join(jobDir,"normalized")
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    cmd = ""
    if method == "UQ":
        outfile = os.path.join(outDir,"matrix_UQ.txt")
        cmd = uq(infile,outfile,jobDir)   
    elif method == "TMM":
        outfile = os.path.join(outDir,"matrix_TMM.txt")
        cmd = tmm(infile,outfile,jobDir)
    elif method == "RLE":
        outfile = os.path.join(outDir,"matrix_RLE.txt")
        cmd = rle(infile,outfile,jobDir)
    elif method == "DESEQ":
        outfile  = os.path.join(outDir,"matrix_DESEQ.txt")
        cmd = deseq(infile,outfile,jobDir)
    elif method == "QN":
        outfile = os.path.join(outDir,"matrix_QN.txt")
        cmd = qn(infile,outfile,jobDir)
    if cmd!="":
        cmds.append(cmd)

    return cmds,outfile
