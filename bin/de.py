import os,sys
import pandas as pd
import numpy as np
import subprocess
from config import R_PATH
from config import R_SCRIPTS_PATH
import itertools
from scipy.stats import ttest_ind
import math

from upsetplot import from_contents,UpSet
from matplotlib import pyplot as plt
from subprocess import Popen

from plotly.offline import plot
import plotly.graph_objects as go
import plotly.express as px

def createGroupFile(annotation_df,jobDir):
    groups = annotation_df["group"].values.tolist()
    diffGroups = list(set(groups))
    combinations = []
    for subset in itertools.combinations(diffGroups, 2):
        combination = subset[0]+"-"+subset[1]
        combinations.append([subset[0],subset[1]])
    output = open(os.path.join(jobDir,"DE","groups.txt"),'a')
    for element in combinations:
        element = "-".join(element)
        output.write(element+"\n")
    output.close()
    return combinations


def de_R(infile,annotation,combinations,method,FDR,min_t,jobDir,log,status):
    commands = []
    for subset in combinations:
        group1 = subset[0]
        group2 = subset[1]

        #edgeR
        output_edgeR = os.path.join(jobDir,"DE","edgeR_"+group1+"_"+group2+".txt")
        cmd_edgeR = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"edgeR_de.R "+infile+" "+method+" "+annotation+" "+FDR+" "+group1+" "+group2+" "+output_edgeR+" 2>"+jobDir+"/Log.txt"
        commands.append(cmd_edgeR)      

        #deseq
        output_deseq = os.path.join(jobDir,"DE","deseq_"+group1+"_"+group2+".txt")
        cmd_deseq = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"deseq_de.R "+infile+" "+annotation+" "+FDR+" "+min_t+" "+group1+" "+group2+" "+output_deseq+" >"+jobDir+"/Log.txt"
        commands.append(cmd_deseq)

        #noiseq
        output_noiseq = os.path.join(jobDir,"DE","noiseq_"+group1+"_"+group2+".txt")
        cmd_noiseq = R_PATH+" --vanilla "+R_SCRIPTS_PATH+"noiseq_de.R "+infile+" "+method+" "+annotation+" "+FDR+" "+min_t+" "+group1+" "+group2+" "+output_noiseq+" >"+jobDir+"/Log.txt"
        commands.append(cmd_noiseq)

    status.write("<p>2. Differential expression analysis has started</p>")
    log.write("2. Differential expression analysis \n")    
    log.write("### EdgeR test DE analysis in progress\n")
    log.write("### DESeq DE analysis in progress\n")
    log.write("### NOISeq DE analysis in progress\n")
    status.flush()

    
    procs = [ Popen(i,shell=True) for i in commands ]
    for p in procs:
        p.wait()
        
    #Read into df
    output={}
    
    for subset in combinations:
        group1 = subset[0]
        group2 = subset[1]

        output_edgeR = os.path.join(jobDir,"DE","edgeR_"+group1+"_"+group2+".txt")
        output_deseq = os.path.join(jobDir,"DE","deseq_"+group1+"_"+group2+".txt")
        output_noiseq = os.path.join(jobDir,"DE","noiseq_"+group1+"_"+group2+".txt")

        try:
            edgeR = pd.read_table(output_edgeR)[["name","logFC","PValue","FDR"]] 
        except:
            edgeR = pd.DataFrame(columns = ["name","logFC","PValue","FDR"])
        
        try:
            deseq = pd.read_table(output_deseq)[["name","logFC","PValue","FDR"]]
        except:
            deseq = pd.DataFrame(columns = ["name","logFC","PValue","FDR"])
        
        try:
            noiseq = pd.read_table(output_noiseq)[["name","logFC","PValue"]]
        except:
            noiseq =  pd.DataFrame(columns=["name","logFC","PValue"])

        if not edgeR.empty:
            log.write("### EdgeR DE analysis finalized\n")
        else:
            log.write("### EdgeR DE analysis finalized with errors\n")
            #status.write("<p>Differential expression analysis with EdgeR finalized with errors</p>")

        if not deseq.empty:
            log.write("### DESeq DE analysis finalized\n")
        else:
            log.write("### DESeq DE analysis finalized with errors\n")
            #status.write("<p>Differential expression analysis with DESeq finalized with errors</p>")

        if not noiseq.empty:
            log.write("### NOISeq DE analysis finalized\n")
        else:
            log.write("### NOISeq DE analysis finalized with errors\n")
            #status.write("<p>Differential expression analysis with NOISeq finalized with errors</p>")
        status.flush()
        output[subset[0]+"-"+subset[1]] = {"edgeR":edgeR,"deseq":deseq,"noiseq":noiseq}


    return output

def ttest(df,combinations,annotation_df,FDR,output_de,jobDir,log,status):

    log.write("### T test DE analysis in progress\n")
    for subset in combinations:

        group1 = subset[0]
        group2 = subset[1]
        output = open(os.path.join(jobDir,"DE","ttest_"+group1+"_"+group2+".txt"),'a')
        cabecera = "name\t"+group1+"_mean\t"+group2+"_mean\tlogFC\tPValue\n"
        output.write(cabecera)
        samples_group1 = annotation_df[annotation_df['group'] ==group1].index.tolist()
        samples_group2 = annotation_df[annotation_df['group'] ==group2].index.tolist()
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
        output_de[subset[0]+"-"+subset[1]]["ttest"] = output_this

    if output_de:
        log.write("### T test DE analysis finalized\n")
    else:
        log.write("### T test DE analysis finalized with errors\n")
        status.write("<p>Differential expression analysis with T test finalized with errors</p>")
    status.flush()
    return output_de

def consensus(df_output,df,annotation_df,jobDir):


    for comparison in df_output:
        sample1 = comparison.split("-")[0]
        sample2 = comparison.split("-")[1]

        #Get expression
        samplesFilterg1 = annotation_df[annotation_df['group'] == sample1]
        samplesGroup1 = list(samplesFilterg1.index)
        mean_g1 = df[samplesGroup1].mean(axis=1).to_frame()

        samplesFilterg2 = annotation_df[annotation_df['group'] == sample2]
        samplesGroup2 = list(samplesFilterg2.index)
        mean_g2 = df[samplesGroup2].mean(axis=1).to_frame()

        expression = pd.merge(mean_g1, mean_g2, left_index =True ,right_index=True)
        expression['log2FC'] = np.where(expression["0_x"]==0,
                        np.abs(np.log2(expression["0_y"])), 
                        np.abs(np.log2(expression["0_y"]/expression["0_x"])))
        expression = expression.round(decimals=2)
        
        try:
            edgeR = df_output[comparison]["edgeR"].name.tolist()
        except:
            edgeR = []

        try:
            deseq = df_output[comparison]["deseq"].name.tolist()
        except:
            deseq=[]

        try:
            noiseq = df_output[comparison]["noiseq"].name.tolist()
        except:
            noiseq = []

        try:
            ttest = df_output[comparison]["ttest"].name.tolist()
        except:
            ttest = []
        methods = from_contents({'edgeR': edgeR, 'DESeq': deseq, 'NOISeq': noiseq,'T-test': ttest })
        upset = UpSet(methods, subset_size='count', show_counts=True,facecolor="#041C34")
        upset.style_subsets(present=["edgeR", "DESeq","NOISeq"],
                    facecolor="#184b7e",
                    label="Consensus 3 methods")
        upset.style_subsets(present=["T-test", "DESeq","NOISeq"],
                    facecolor="#184b7e",
                    label="Consensus 3 methods")
        upset.style_subsets(present=["edgeR", "T-test","NOISeq"],
                    facecolor="#184b7e",
                    label="Consensus 3 methods")
        upset.style_subsets(present=["edgeR", "T-test","DESeq"],
                    facecolor="#184b7e",
                    label="Consensus 3 methods")
        upset.plot()
        plt.savefig(os.path.join(jobDir,"DE","consensus_upset_"+sample1+"-"+sample2+".png"))


        edgeR = df_output[comparison]["edgeR"]
        deseq = df_output[comparison]["deseq"]
        noiseq = df_output[comparison]["noiseq"]
        ttest = df_output[comparison]["ttest"]

        try:
            result = pd.merge(edgeR, deseq, on="name",how="outer")
            result.rename(columns = {'FDR_x':'edgeR','FDR_y':'DESeq'}, inplace = True)
            result = result[['name','edgeR','DESeq']]
            result = pd.merge(result, noiseq, on="name",how="outer", suffixes= ["","_noiseq"])
            result.rename(columns = {'PValue':'NOISeq'}, inplace = True)
            result = result[['name','edgeR','DESeq','NOISeq']]
            result = pd.merge(result, ttest, on="name",how="outer", suffixes= ["","_ttest"])
            result.rename(columns = {'PValue':'T-Test',}, inplace = True)
            result = result[['name','edgeR','DESeq','NOISeq','T-Test']]
            result["edgeR"].fillna("Non-DE",inplace=True)
            result["DESeq"].fillna("Non-DE",inplace=True)
            result["NOISeq"].fillna("Non-DE",inplace=True)
            result["T-Test"].fillna("Non-DE",inplace=True)
            result = result.set_index('name')
            result_filtered = pd.merge(result, expression, left_index =True,right_index=True)
            result_filtered.rename(columns = {'0_x':sample1+"_mean",'0_y':sample2+"_mean"}, inplace=True)
            numMethods = 4 - (result_filtered.isin(['Non-DE']).sum(axis=1)).to_frame()
            dfConsensus = pd.merge(result_filtered, numMethods, left_index =True,right_index=True)
            dfConsensus.rename(columns = {0:"Number Methods"}, inplace=True)
            outfile_consensus = os.path.join(jobDir,"DE","consensus_"+sample1+"-"+sample2+".txt") 
            dfConsensus.to_csv(outfile_consensus,sep="\t")
        except:
            pass

def plotDE(df,output_de,annotation_df,jobDir):
    for combination in output_de:
        for method in output_de[combination]:
            outfile = os.path.join(jobDir,"DE","top10_"+method+"_"+combination+".html") 
            myDf = output_de[combination][method]
            myDf = myDf.set_index("name")
            try:
                top10 = list((myDf.nlargest(10, 'logFC')).index)
                result_expr = df.loc[top10].T
                result_expr = pd.merge(result_expr, annotation_df, left_index =True ,right_index=True)
                fig = go.Figure()

                fig = px.box(result_expr,color='group',
                labels={
                                "group": "Group",
                                "value":"Expression",
                                "variable":""
                            },)

                
                plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')
                outfile_W = open(outfile,'a')
                outfile_W.write(plotCode)
                outfile_W.close()
            except:
                pass
