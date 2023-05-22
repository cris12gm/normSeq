import os,sys
import pandas as pd
import numpy as np
import subprocess

from config import R_PATH
from config import R_SCRIPTS_PATH

def processInputInit(infile,samples,minRC,annotation_df):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df=df.dropna(axis=1,how='all')
    test = df.select_dtypes(include=[float])

    if test.empty:
        df = pd.read_table(infile,decimal=",")
        df=df.dropna(axis=1,how='all')

    secondtest = df.select_dtypes(include=[float])
    if secondtest.empty:
        sys.exit(0)

    df = df.rename(columns=lambda s: s.replace("_", ""))
    df = df.rename(columns=lambda s: s.replace(" ", ""))

    if cabecera!="":
        df.rename(columns = {cabecera:'name'}, inplace = True)
    else:
        df.rename(columns = {list(df)[0]:'name'}, inplace=True)
    df = df.set_index('name')
    
    samples_df = list(set(samples) & set(df.columns))
    df = df.dropna()
    dfF = df[samples_df]
    try:
        dfF = dfF[dfF>=0]
    except:
        df = pd.read_table(infile,decimal=",")
        df = df.rename(columns=lambda s: s.replace("_", ""))
        df = df.rename(columns=lambda s: s.replace(" ", ""))
        if cabecera!="":
            df.rename(columns = {cabecera:'name'}, inplace = True)
        else:
            df.rename(columns = {list(df)[0]:'name'}, inplace=True)
        df = df.set_index('name')
        df = df.dropna()
        dfF = df[samples]

    dfOriginal = dfF
    minRC = int(minRC)
    dfF = dfF.loc[~(dfF==0).all(axis=1)]
    dfF = dfF[dfF>=minRC]
    dfF=dfF.dropna(axis=0)

    try:
        samples = annotation_df.index.values
        dfF=dfF[samples]
    except:
        pass

    return(dfF,dfOriginal)


def processInput(infile,annotation_df):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df = df.rename(columns=lambda s: s.replace("_", ""))
    df = df.rename(columns=lambda s: s.replace(" ", ""))
    if cabecera!="":
        df.rename(columns = {cabecera:'name'}, inplace = True)
    else:
        df.rename(columns = {list(df)[0]:'name'}, inplace=True)
    df = df.set_index('name')
    dfF = df.dropna()

    samples = annotation_df.index.values
    dfF=dfF[samples]
    
    return(dfF)

def processAnnotation(infile):
    cabecera = open(infile).readline().split("\t")
    df = pd.read_table(infile,index_col=False)
    numbercolumns = len(df.columns)
    df = df.replace(' ','', regex=True)
    df = df.replace('_','', regex=True)
    
    
    if numbercolumns == 3:
        df.rename(columns = {cabecera[0]:'sample',cabecera[1]:'group',cabecera[2]:'replicate'}, inplace = True)
        df.columns = ['sample','group','replicate']
        replicates = True
    elif numbercolumns == 2:
        df.rename(columns = {cabecera[0]:'sample',cabecera[1]:'group'}, inplace = True)
        df.columns = ['sample','group']
        replicates = False
    else:
        df = df.iloc[:, 0:2]
        df.rename(columns = {cabecera[0]:'sample',cabecera[1]:'group'}, inplace = True)
        df.columns = ['sample','group']
        replicates = False
    df = df.set_index('sample')
    df.dropna(how='all', axis=1, inplace=True)
    samples = df.index.tolist()
    return(df,samples,replicates)

def cpm(df,outfile):
    df = df.loc[~(df==0).all(axis=1)]
    #Normalize CPM
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    normalized = np.divide(df, sums)*1000000

    normalized.to_csv(outfile,sep="\t")

    return normalized

# def tc(df,outfile):
#     sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
#     mean = np.mean(sums)

#     normalized = np.divide(df,mean)*1000000
#     normalized.to_csv(outfile,sep="\t")

#     return normalized

def uq(df,outfile):
    df = df.loc[~(df==0).all(axis=1)]
    #Execute in R
    # Calculate the upper quartile normalization factors
    upper_quartile = np.percentile(df, 75, axis=1)

    # Replace any zero normalization factors with median non-zero value
    med_norm_factor = np.median(upper_quartile[upper_quartile > 0])
    upper_quartile[upper_quartile == 0] = med_norm_factor

    # Divide the counts by the normalization factors
    norm_factors = df.divide(upper_quartile, axis=0)

    # Save the normalized counts to a new file
    normalized = norm_factors.multiply(np.median(upper_quartile))
    normalized.to_csv(outfile,sep="\t")

    return normalized

def tmm(infile,outfile,jobDir):
    #Execute in R
    cmd_tmm = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_normalization.R "+infile+" TMM "+outfile+" >"+jobDir+"/Log.txt 2>"+jobDir+"/Log.txt"
    return cmd_tmm

def rle(infile,outfile,jobDir):
    #Execute in R
    cmd_rle = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_normalization.R "+infile+" RLE "+outfile+" >"+jobDir+"/Log.txt 2>"+jobDir+"/Log.txt"
    return cmd_rle

def med(df,outfile):
    df = df.loc[~(df==0).all(axis=1)]
    sums = np.array(pd.DataFrame(np.sum(df,axis=0)).T)
    median = np.median(sums)

    normalized = np.divide(df,median)*1000000
    normalized.to_csv(outfile,sep="\t")

    return normalized

def deseq(infile,outfile,jobDir):
    #Execute in R
    cmd_deseq = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"deseq_normalization.R "+infile+" "+outfile+" >"+jobDir+"/Log.txt 2>"+jobDir+"/Log.txt"
    return cmd_deseq


def qn(infile,outfile,jobDir):
    #Execute in R
    cmd_qn = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"quantile_normalization.R "+infile+" "+outfile+" >"+jobDir+"/Log.txt 2>"+jobDir+"/Log.txt"
    return cmd_qn

def ruv(infile,outfile,annotation,jobDir):
    cmd_ruv = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"ruv_normalization.R "+infile+" "+annotation+" "+outfile+" >"+jobDir+"/Log.txt 2>"+jobDir+"/Log.txt"
    return cmd_ruv

def norm(infile,df,method,jobDir):
    outDir = os.path.join(jobDir,"normalized")
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    if method == "CPM":
        outfile = os.path.join(outDir,"matrix_CPM.txt")
        outdf = cpm(df,outfile)
    # elif method =="TC":
    #     outfile = os.path.join(outDir,"matrix_TC.txt")
    #     outdf = tc(df,outfile)
    elif method == 'NN':
        outfile = os.path.join(outDir,"matrix_NN.txt") 
        outdf = df
    elif method == 'Med':
        outfile = os.path.join(outDir,"matrix_Med.txt")
        outdf = med(df,outfile)
    elif method == 'UQ':
        outfile = os.path.join(outDir,"matrix_UQ.txt")
        outdf = uq(df,outfile)

    return outdf,outfile


def norm_r(infile,method,jobDir,annotation,cmds):
    outDir = os.path.join(jobDir,"normalized")
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    cmd = ""

    if method == "TMM":
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
    elif method == "RUV":
        outfile = os.path.join(outDir,"matrix_RUV.txt")
        cmd = ruv(infile,outfile,annotation,jobDir)
    if cmd!="":
        cmds.append(cmd)

    return cmds,outfile


def rleplot(infile,annotation,jobDir,method):
    outfile = os.path.join(jobDir,"normalized","RLEplot_"+method+".png")
    cmd_rle = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"plotRLE.R "+infile+" "+annotation+" "+outfile+" >"+jobDir+"/Log.txt"
    return cmd_rle