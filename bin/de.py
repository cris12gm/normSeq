import os,sys
import pandas as pd
import numpy as np
import subprocess
from config import R_PATH
from config import R_SCRIPTS_PATH
import itertools
from scipy.stats import ttest_ind
import math

import upsetplot
from matplotlib import pyplot

def processAnnotation(infile):
    cabecera = open(infile).readline().split("\t")[0]
    df = pd.read_table(infile)
    df.rename(columns = {cabecera:'sample'}, inplace = True)
    df = df.set_index(cabecera)
    df = df.dropna()
    return(df)

def createGroupFile(annotation,jobDir):
    annotationFile = processAnnotation(annotation)
    groups = annotationFile["group"].values.tolist()
    diffGroups = list(set(groups))
    combinations = []
    for subset in itertools.combinations(diffGroups, 2):
        combination = subset[0]+"-"+subset[1]
        combinations.append(combination)
    output = open(os.path.join(jobDir,"DE","groups.txt"),'a')
    for element in combinations:
        output.write(element+"\n")
    output.close()
    return combinations


def edgeR(infile,method,annotation,FDR,jobDir):

    annotationFile = processAnnotation(annotation)
    groups = annotationFile["group"].values.tolist()
    diffGroups = list(set(groups))
    outputdf = {}

    for subset in itertools.combinations(diffGroups, 2):
        group1 = subset[0]
        group2 = subset[1]
        output = os.path.join(jobDir,"DE","edgeR_"+group1+"_"+group2+".txt")
        subprocess.call (R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_de.R "+infile+" "+method+" "+annotation+" "+FDR+" "+group1+" "+group2+" "+output,shell=True)

        output_this = pd.read_table(output)[["name","logFC","PValue","FDR"]]
        outputdf[group1,group2] = output_this

    return outputdf

def deseq(infile,annotation,FDR,min_t,jobDir):

    annotationFile = processAnnotation(annotation)
    groups = annotationFile["group"].values.tolist()
    diffGroups = list(set(groups))

    outputdf = {}
    for subset in itertools.combinations(diffGroups, 2):
        group1 = subset[0]
        group2 = subset[1]
        output = os.path.join(jobDir,"DE","deseq_"+group1+"_"+group2+".txt")
        subprocess.call (R_PATH+" --vanilla "+R_SCRIPTS_PATH+"deseq_de.R "+infile+" "+annotation+" "+FDR+" "+min_t+" "+group1+" "+group2+" "+output,shell=True)

        #FIX!
        try:
            output_this = pd.read_table(output)[["name","logFC","PValue","FDR"]]
            outputdf[group1,group2] = output_this
        except:
            pass
    return outputdf

def noiseq(infile,method,annotation,FDR,min_t,jobDir):

    annotationFile = processAnnotation(annotation)
    groups = annotationFile["group"].values.tolist()
    diffGroups = list(set(groups))
    outputdf = {}
    for subset in itertools.combinations(diffGroups, 2):
        group1 = subset[0]
        group2 = subset[1]
        output = os.path.join(jobDir,"DE","noiseq_"+group1+"_"+group2+".txt")
        subprocess.call (R_PATH+" --vanilla "+R_SCRIPTS_PATH+"noiseq_de.R "+infile+" "+method+" "+annotation+" "+FDR+" "+min_t+" "+group1+" "+group2+" "+output,shell=True)

        output_this = pd.read_table(output)[["name","logFC","PValue"]]
        outputdf[group1,group2] = output_this

    return outputdf

def ttest(df,annotation,FDR,jobDir):
    annotationFile = processAnnotation(annotation)
    groups = annotationFile["group"].values.tolist()
    diffGroups = list(set(groups))

    outputdf = {}
    for subset in itertools.combinations(diffGroups, 2):

        group1 = subset[0]
        group2 = subset[1]
        output = open(os.path.join(jobDir,"DE","ttest_"+group1+"_"+group2+".txt"),'a')
        cabecera = "name\t"+group1+"_mean\t"+group2+"_mean\tlogFC\tPValue\n"
        output.write(cabecera)
        samples_group1 = annotationFile[annotationFile['group'] ==group1].index.tolist()
        samples_group2 = annotationFile[annotationFile['group'] ==group2].index.tolist()
        output_this = {'name':[],'logFC':[],'PValue':[]}


        for mirna, expression in df.iterrows():
            group1_values = expression[samples_group1].tolist()
            group2_values = expression[samples_group2].tolist()
            ttest = ttest_ind(group1_values, group2_values,
                      equal_var=False, # it's not necessarily fair to assume that these two populations have equal variance
                      nan_policy='omit') # omit NaN values
            if float(ttest.pvalue)<=float(FDR):
                mean_group1 = round(sum(group1_values) / len(group1_values),2)
                mean_group2 = round(sum(group2_values) / len(group2_values),2)
                log2 = math.log((mean_group1+1)/(mean_group2+1),2)
                ttest = str(round(ttest.pvalue,2))
                escribir = mirna+"\t"+str(mean_group1)+"\t"+str(mean_group2)+"\t"+str(log2)+"\t"+ttest+"\n"
                output.write(escribir)
                output_this['name'].append(mirna)
                output_this['logFC'].append(log2)
                output_this['PValue'].append(ttest)
                output_this[group1+"_mean"] = mean_group1
                output_this[group2+"_mean"] = mean_group2
        output.close()
        output_this = pd.DataFrame(output_this)
        outputdf[group1,group2] = output_this
    return outputdf
def consensus(df1,df2,df3,df4,jobDir):
    out = {}
    for comparison in df4:
        edgeR = df1[comparison]
        deseq = df2[comparison]
        noiseq = df3[comparison]
        ttest = df4[comparison]
        result = pd.merge(edgeR, deseq, on="name",how="outer")
        result = result[['logFC_x','logFC_y']]
        result.rename(columns = {'logFC_x':'edgeR','logFC_y':'deseq'}, inplace = True)
        result["edgeR"].fillna(False,inplace=True)
        result["deseq"].fillna(False,inplace=True)
        result.edgeR = result.edgeR.astype('bool')
        result.deseq = result.deseq.astype('bool')
        plotUpset(result)
        print (result)
        sys.exit(1)
        result = pd.merge(result, ttest, on="name",how="outer", suffixes= ["","_ttest"])
        result = pd.merge(result, noiseq, on="name",how="outer", suffixes= ["","_noiseq"])
        result_filtered = result[['PValue_edgeR', 'PValue_deseq',"PValue","PValue_noiseq"]]
        result_filtered.rename(columns = {'PValue_edgeR':'edgeR','PValue_deseq':'deseq','PValue':'ttest','PValue_noiseq':'noiseq'}, inplace = True)
        result_filtered.edgeR_2 = result_filtered.edgeR.astype('bool')
        print (result_filtered)
#        out[comparison] = result

#        plotUpset(result)

