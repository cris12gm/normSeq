from django.conf import settings
import random,string
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic.edit import CreateView
from django.views.generic import FormView, DetailView, TemplateView

from django.http import JsonResponse

import os
import pandas as pd
import json


METHODS = {
        'NN' :' No normalization',
        'CPM' : ' Counts per Million',
        'TC' : ' Total Count',
        'UQ' :'Upper Quartile',
        'Med' :' Median',
        'DESEQ' : 'DEseq',
        'TMM' : ' TMM',
        'QN' : ' Quantile',
        'RUV' : ' Remove Unwanted Variation',
        'RLE' : 'Relative Log Expression'
    }
class Errors(Enum):
    NO_ERROR = 0
    NOT_VALID = 1
    NOT_ASSOCIATED = 2

class Results(TemplateView):
    template = 'results.html'

    def get(self, request):

        if not request.GET['jobID']:
            return redirect(settings.SUB_SITE+"/query/")

        jobID = request.GET['jobID']

        #Get config file

        configFile = open(settings.MEDIA_ROOT+jobID+'/config.json')
        config = json.load(configFile)

        methods = config['methods']
        jobDir = config['jobDir']
        typeJob = config['typeJob']

        ##All visualizations
        visualization = False

        heatmap = []
        pca = []
        downloads = {}
        distribution = {}
        top10 = {}
        de_groups={}
        features = {}

        i = 0

        #InfoGain
        infoGain = {}
        info = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain.html")
        infoPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain.png")
        infoGain['HTML'] = info
        infoGain['PNG'] = infoPNG
        infoGain['name'] = "Information Gain"

        ##Dif expr
        diffExpr = {}
        if config['diffExpr'] == 'True':
            groupsFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"DE","groups.txt"),'r')
            groups = []
            for line in groupsFile:
                groups.append(line.strip())
            de_groups = groups
                
            for group in groups:
                group = group.replace("-","_")
                group1 = group.split("_")[0]
                group2 = group.split("_")[1]
                summaryDE,resultsGroup = de_prepare(jobID,group,group1,group2)
                group = group.replace("_","-")
                diffExpr[group]= summaryDE,resultsGroup
                
        #Visualization
        for method in methods:

            pngHeatmap = os.path.join(settings.MEDIA_URL,jobID,"graphs","heatmap_"+method+".png")
            id_modal = "heatmap_"+method
            title_modal = METHODS[method]
            #heatmapHTML = os.path.join(settings.MEDIA_URL,jobID,"graphs","heatmap_"+method+".html")

            # heatmap.append([pngHeatmap,heatmapHTML,id_modal,title_modal])
            heatmap.append([pngHeatmap,id_modal,title_modal])
            

            pngPCA = os.path.join(settings.MEDIA_URL,jobID,"graphs","pca_"+method+".png")
            id_modal = "pca_"+method
            title_modal = METHODS[method]
            pcaHTML = os.path.join(settings.MEDIA_URL,jobID,"graphs","pca_"+method+".html")

            pca.append([pngPCA,pcaHTML,id_modal,title_modal])

            #Summary
            #Distribution
            distribution[method] = {}
            distributionHTML = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","distribution_"+method+".html")
            distributionPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","distribution_"+method+".png")
            distribution[method]['HTML'] = distributionHTML
            distribution[method]['PNG'] = distributionPNG
            distribution[method]['name'] = METHODS[method]
            if i==0:
                distribution[method]['active'] = "block;"
            else:
                distribution[method]['active'] = "none;"
               

            #Top miRNAS
            top10[method] = {}
            top = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10_"+method+".html")
            topPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10_"+method+".png")
            top10[method]['HTML'] = top
            top10[method]['PNG'] = topPNG
            top10[method]['name'] = METHODS[method]
            if i==0:
                top10[method]['active'] = "block;"
            else:
                top10[method]['active'] = "none;"

            if config['batchEffect'] == 'True':
                #Batch effect Plots
                batchEffect = {}
                before = os.path.join(settings.MEDIA_URL,jobID,"graphs","batchEffect","pca_old.html")
                beforePNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","batchEffect","pca_old.png")
                after = os.path.join(settings.MEDIA_URL,jobID,"graphs","batchEffect","pca_corrected.html")
                afterPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","batchEffect","pca_corrected.png")
            
                batchEffect['before'] = before
                batchEffect['beforePNG'] = beforePNG
                batchEffect['after'] = after
                batchEffect['afterPNG'] = afterPNG
            else:
                batchEffect=False


            ##Downloads
            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"normalized","matrix_"+method+".txt")

            downloads[METHODS[method]] = [downloadLink,"matrix_"+method+".txt"]
            i = i + 1 

        summary = {'distribution':distribution,'top10':top10,'info':infoGain}

        if heatmap or pca:
            visualization=True

        return render(request, self.template, {"jobID":jobID,"typeJob":typeJob,"visualization":visualization,"heatmapPlots":heatmap,
        "pcaPlots":pca,"batchEffect":batchEffect,"downloads":downloads,"summary":summary,"de":diffExpr,"de_groups":de_groups})


def queryPlotHTML(request):
#    templateError = "error.html"
    # try:
    url = request.GET.get('url', None).replace(settings.MEDIA_URL,settings.MEDIA_ROOT)
    with open(url,'r') as file:
        plot = file.read().rstrip()

    data = {}
    data["plot"]=plot

    return JsonResponse(data)
    # except:
    #     return render(request, templateError)


def queryPlotFeature(request):
    method = request.GET.get('method', None)
    print(method)
    # url = request.GET.get('url', None).replace(settings.MEDIA_URL,settings.MEDIA_ROOT)
    # with open(url,'r') as file:
    #     plot = file.read().rstrip()

    data = {}
    # data["plot"]=plot

    return JsonResponse(data)
    # except:
    #     return render(request, templateError)

def de_prepare(jobID,group,group1,group2):
    summary = {}
    resultsGroup = {}
    #edgeR
    edgeRFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","edgeR_"+group+".txt")
    try:
        edgeRdf = pd.read_table(edgeRFile)
        edgeR = edgeRdf[["name",group1+"_mean", group2+"_mean","logFC","PValue","FDR"]]
        edgeR = edgeR.set_index("name")
        edgeR = edgeR.round(decimals=3)
        summary['edgeR'] = getInfraOver(edgeR)

        edgeR_header = list(edgeR.columns)
        edgeR = edgeR.T.to_dict()
        resultsGroup['edgeR_header'] = edgeR_header
        resultsGroup['edgeR'] = edgeR
    except:
        pass
        
    #deseq
    deseqFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","deseq_"+group+".txt")
    try:
        deseqdf = pd.read_table(deseqFile)
        deseq = deseqdf[["name",group1+"_mean", group2+"_mean","logFC","PValue","FDR"]]
        deseq = deseq.set_index("name")
        deseq = deseq.round(decimals=3)
        summary['DESeq2'] = getInfraOver(deseq)
        deseq_header = list(deseq.columns)
        deseq = deseq.T.to_dict()
        resultsGroup['deseq_header'] = deseq_header
        resultsGroup['deseq'] = deseq
    except:
        pass
    
    #noiseq
    noiseqFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","noiseq_"+group+".txt")
    try:
        noiseqdf = pd.read_table(noiseqFile)
        noiseq = noiseqdf[["name",group1+"_mean", group2+"_mean","logFC","PValue"]]
        noiseq = noiseq.set_index("name")
        noiseq = noiseq.round(decimals=3)
        summary['NOISeq'] = getInfraOver(noiseq)
        noiseq_header = list(noiseq.columns)
        noiseq = noiseq.T.to_dict()
        resultsGroup['noiseq_header'] = noiseq_header
        resultsGroup['noiseq'] = noiseq
    except:
        pass
    #ttest
    ttestFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","ttest_"+group+".txt")
    try:
        ttestdf = pd.read_table(ttestFile)
        ttest = ttestdf[["name",group1+"_mean", group2+"_mean","logFC","PValue"]]
        ttest = ttest.set_index("name")
        summary['T Test'] = getInfraOver(ttest)
        ttest_header = list(ttest.columns)
        ttest = ttest.T.to_dict()
        resultsGroup['ttest_header'] = ttest_header
        resultsGroup['ttest'] = ttest
    except:
        pass
    return summary,resultsGroup

def getInfraOver(df):
    output = {}

    overExpressed = df.query("logFC >= 0")
    overExpressed_number = overExpressed.shape[0]
    overExpressed_names = overExpressed.index.tolist()
    overExpressed_names = ','.join(overExpressed_names)

    infraExpressed = df.query("logFC < 0")
    infraExpressed_number = infraExpressed.shape[0]
    infraExpressed_names = infraExpressed.index.tolist()
    infraExpressed_names = ','.join(infraExpressed_names)

    output = {"over":overExpressed_names,"infra":infraExpressed_names,"numOver":overExpressed_number,"numInfra":infraExpressed_number}

    return output
