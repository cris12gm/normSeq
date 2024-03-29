from django.conf import settings
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic.edit import CreateView
from django.views.generic import FormView, DetailView, TemplateView

from django.http import JsonResponse

import os,sys
import pandas as pd
import json
import datetime
import plotly.express as px
from plotly.offline import plot
import plotly.graph_objects as go
from query.models import Job


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

        jobDB = Job.objects.get(JobID=jobID)
        typeJob = jobDB.getType()
        start_date = jobDB.getDate()
        end_date = start_date + datetime.timedelta(days=15)


        #Get warning file

        
        if os.path.exists (settings.MEDIA_ROOT+jobID+'/warnings.txt'):
            warningFileContent = open(settings.MEDIA_ROOT+jobID+'/warnings.txt','r').readlines()[0]
        else:
            warningFileContent=False
        ##All visualizations
        visualization = False

        heatmap = []
        pca = []
        downloads = {}
        distribution = {}
        top10 = {}
        top10FC = {}
        de_groups={}
        features = {}
        rleplots = {}

        i = 0
        j = 0

        Info_Groups = []
        try:
            groupsFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","summary","groups.txt"),'r')
            for line in groupsFile:
                Info_Groups.append(line.strip())

        except:
            pass

        FC_groups = []
        try:
            groupsFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","summary","combinations.txt"),'r')
            for line in groupsFile:
                FC_groups.append(line.strip())
        except:
            pass


        #InfoGain per group
        infoGain = {}
        downloads['info'] = {}
        for group in Info_Groups:
            info = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".html")
            if os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","summary","infoGain_"+group+".html")):
                infoPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".png")
                infoGain[group] = {'HTML':info,'PNG':infoPNG,'name':"Information Gain "}
                downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".txt")
                try:
                    values = downloads['info'][group]
                except:
                    downloads['info'][group] = {}
                downloads['info'][group] = [downloadLink,"infoGain_"+group+".txt"]
        infoGainPairwise = {}
        for group in FC_groups:
            info = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".html")
            if os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","summary","infoGain_"+group+".html")):
                infoPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".png")
                infoGainPairwise[group] = {'HTML':info,'PNG':infoPNG,'name':"Information Gain "}

                downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".txt")
                try:
                    values = downloads['info'][group]
                except:
                    downloads['info'][group] = {}
                downloads['info'][group] = [downloadLink,"infoGain_"+group+".txt"]
        

        ## Features
        no_normalized = os.path.join(settings.MEDIA_ROOT,jobID,"normalized","matrix_NN.txt")
        features = (pd.read_table(no_normalized)['name']).tolist()
        features.sort()

        ##Dif expr
        diffExpr = {}
        consensus = {}
        topDEPerMethod = {}
        downloads['DE']={}
        if config['diffExpr'] == 'True':
            try:
                groupsFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"DE","groups.txt"),'r')
                groups = []
                for line in groupsFile:
                    groups.append(line.strip())
                de_groups = groups
                
                for group in groups:
                    group = group.replace("-","_")
                    groupFile = group
                    group1 = group.split("_")[0]
                    group2 = group.split("_")[1]
                    summaryDE,resultsGroup,consensusM = de_prepare(jobID,group,group1,group2)
                    group = group.replace("_","-")
                    diffExpr[group]= summaryDE,resultsGroup
                    consensus[group] = consensusM
                    topDEPerMethod[group] = de_prepareTop(jobID,group)

                    #Downloads DE
                    protocols = {}
                    protocols = {"edgeR":"edgeR","DESeq2":"deseq2","NOISeq":"noiseq","T-Test":"test"}
                    
                    for protocolName in protocols:
                        protocolFile = protocols[protocolName]
                        downloadLink = os.path.join(settings.MEDIA_URL,jobID,"DE",protocolFile+"_"+groupFile+".txt")
                        try:
                            values = downloads['DE'][protocolName]
                        except:
                            downloads['DE'][protocolName] = {}
                        downloads['DE'][protocolName][group] = [downloadLink,protocolFile+"_"+groupFile+".txt"]
            except:
                pass

        ##Downloads All
        downloadLink = os.path.join(settings.MEDIA_URL,jobID,"results_"+jobID+".zip")
        downloads['all'] = {}
        downloads['all'] = [downloadLink,"results_"+jobID+".zip"]

        downloads['normalized'] = {}
        downloads['top10FC'] = {}
        downloads['top10'] = {}
        #Visualization
        for method in methods:

            pngHeatmap = os.path.join(settings.MEDIA_URL,jobID,"graphs","heatmap_"+method+".png")
            if os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","heatmap_"+method+".png")):
                id_modal = "heatmap_"+method
                title_modal = METHODS[method]
                #heatmapHTML = os.path.join(settings.MEDIA_URL,jobID,"graphs","heatmap_"+method+".html")

                # heatmap.append([pngHeatmap,heatmapHTML,id_modal,title_modal])
                heatmap.append([pngHeatmap,id_modal,title_modal])
            

            pngPCA = os.path.join(settings.MEDIA_URL,jobID,"graphs","pca_"+method+".png")
            pngPCA3D = True
            if os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","pca_"+method+".png")):
                id_modal = "pca_"+method
                title_modal = METHODS[method]
                pcaHTML = os.path.join(settings.MEDIA_URL,jobID,"graphs","pca_"+method+".html")

                pca.append([pngPCA,pcaHTML,id_modal,title_modal])

                #Check if 3D plot
            pngPCA3D = os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","pca_"+method+"_3D.html"))

            #Summary
            #Distribution
            if os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","summary","distribution_"+method+".html")):
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

            #Downloads
            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10_"+method+".txt")
            try:
                values = downloads['top10'][METHODS[method]]
            except:
                downloads['top10'][METHODS[method]] = []
            downloads['top10'][METHODS[method]] = [downloadLink,method+".txt"]

            #Top miRNAS FC
            top10FC[method] = {}

            if FC_groups:

                top10FC[method]['name'] = METHODS[method]
                top10FC[method]["comparisons"] = FC_groups
                top10FC[method]["data"] = {}

                for element in FC_groups:
                    topFC = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10FC_"+method+"_"+element+".html")
                    topFCPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10FC_"+method+"_"+element+".png")
                    try:
                        value = top10FC[method]["data"][element]
                    except:
                        top10FC[method]["data"][element] = {}
                    top10FC[method]["data"][element]['HTML'] = topFC
                    top10FC[method]["data"][element]['PNG'] = topFCPNG
                    top10FC[method]["data"][element]['name'] = METHODS[method]
                    if j==0:
                        top10FC[method]["data"][element]['active'] = "block;"
                    else:
                        top10FC[method]["data"][element]['active'] = "none;"

                    #Downloads
                    downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10FC_"+method+"_"+element+".txt")
                    try:
                        values = downloads['top10FC'][METHODS[method]]
                    except:
                        downloads['top10FC'][METHODS[method]] = {}
                    downloads['top10FC'][METHODS[method]][element] = [downloadLink,method+"_"+element+".txt"]

                    j = 1
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

            #RLE plots
            if method=="NN":
                NN_PNG = os.path.join(settings.MEDIA_URL,jobID,"normalized","RLEplot_NN.png")
                rleplotNN={}
                if os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"normalized","RLEplot_NN.png")):
                    rleplotNN['PNG'] = NN_PNG
                    rleplotNN['name'] = METHODS[method]
            else:
                rlePNG = os.path.join(settings.MEDIA_URL,jobID,"normalized","RLEplot_"+method+".png")
                if os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID,"normalized","RLEplot_"+method+".png")):
                    rleplots[method] = {}
                    rleplots[method]['PNG'] = rlePNG
                    rleplots[method]['name'] = METHODS[method]
                    if i==1:
                        rleplots[method]['active'] = "block;"
                    else:
                        rleplots[method]['active'] = "none;"
            try:
                value = rleplotNN
            except:
                rleplotNN={}
            ##Downloads Normalized
            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"normalized","matrix_"+method+".txt")
            downloads['normalized'][METHODS[method]] = [downloadLink,"matrix_"+method+".txt"]

            i = i + 1 

        summary = {'distribution':distribution,'top10':top10,'top10FC':top10FC,'FC_groups':FC_groups,'info':infoGain,'infoPairwise':infoGainPairwise,'rleplotNN':rleplotNN,'rleplots':rleplots}

        if heatmap or pca:
            visualization=True

        de_software = ["edgeR","DESeq2","NOISeq","TTest"]
        return render(request, self.template, {"jobID":jobID,"typeJob":typeJob,"visualization":visualization,"heatmapPlots":heatmap,
        "pcaPlots":pca,"pngPCA3D":pngPCA3D,"batchEffect":batchEffect,"downloads":downloads,"summary":summary,"de":diffExpr,"de_groups":de_groups,'date':end_date,
        'consensus':consensus,'de_software':de_software, 'topDEPerMethod':topDEPerMethod,'features':features,'warning':warningFileContent})


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

def queryPlotTopDE(request):
    method = request.GET.get('method')
    comparison = request.GET.get('comparison')
    jobID = request.GET.get('jobID')
    
    if method=="T Test":
        method="ttest"
    elif method=="DESeq2":
        method ="deseq"

    url = os.path.join(settings.MEDIA_ROOT,jobID,"DE","top10_"+method+"_"+comparison+".html")
    with open(url,'r') as file:
        plot = file.read().rstrip()
    data = {}
    data["plot"]=plot
    return JsonResponse(data)



def queryPlotFeature(request):
    feature = request.GET.get('feature', None)
    jobID = request.GET.get('jobID', None)
    configFile = open(settings.MEDIA_ROOT+jobID+'/config.json')
    config = json.load(configFile)

    methods = config['methods']
    
    annotationdf = pd.read_table(os.path.join(settings.MEDIA_ROOT,jobID,"annotation.txt"))
    annotationdf = annotationdf.set_index('sample')
    toPlot = []
    for method in methods:
        normDf = pd.read_table(os.path.join(settings.MEDIA_ROOT,jobID,"normalized","matrix_"+method+".txt"))
        normDf = normDf.set_index('name')

        thisSample = normDf.loc[feature]
        thisSample = pd.merge(thisSample,annotationdf,right_index=True,left_index=True)
        thisSample['method'] = METHODS[method]
        toPlot.append(thisSample)
    result = pd.concat(toPlot)

    fig = px.box(result, x="group", y=feature, color="group",facet_col="method", facet_col_wrap=4,facet_col_spacing=0.04)
    fig.update_yaxes(matches=None,showticklabels=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_layout(
    width=1400,
    height=500,
    font=dict(
        size=16,
    )
    )
    for axis in fig.layout:
        if type(fig.layout[axis]) == go.layout.XAxis:
            fig.layout[axis].title.text = ''


    plotCode = plot(fig, show_link=False, auto_open=False, output_type = 'div')

    data = {}
    data["plot"]=plotCode

    return JsonResponse(data)

def de_prepareTop(jobID,group):
    resultsGroup =  {}
    #edgeR
    edgeRFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","top10_edgeR_"+group+".html")
    resultsGroup["edgeR"] = edgeRFile
    #DESEq2
    deseqFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","top10_deseq_"+group+".html")
    resultsGroup["DESeq2"] = deseqFile
    #NOISeq
    noiseqFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","top10_noiseq_"+group+".html")
    resultsGroup["NOISeq"] = noiseqFile
    #T Test
    ttestFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","top10_ttest_"+group+".html")
    resultsGroup["TTest"] = ttestFile

    return(resultsGroup)

    

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

    ##Consensus
    consensus={}
    groupC = group.replace("_","-")
    consensusFile = os.path.join(settings.MEDIA_ROOT,jobID,"DE","consensus_"+groupC+".txt")
    try:
        consensusdf = pd.read_table(consensusFile)
        consensusdf = consensusdf.set_index("name")
        consensus_header = list(consensusdf.columns)
        consensusTable = consensusdf.T.to_dict()
        consensus['header'] = consensus_header
        consensus['table'] = consensusTable
        plot = os.path.join(settings.MEDIA_URL,jobID,"DE","consensus_upset_"+groupC+".png")
        consensus['plot'] = plot
    except:
        pass
    return summary,resultsGroup,consensus

def getInfraOver(df):
    output = {}

    overExpressed = df.query("logFC <= 0")
    overExpressed_number = overExpressed.shape[0]
    overExpressed_names = overExpressed.index.tolist()
    overExpressed_names = ','.join(str(e) for e in overExpressed_names)
    infraExpressed = df.query("logFC > 0")
    infraExpressed_number = infraExpressed.shape[0]
    infraExpressed_names = infraExpressed.index.tolist()
    infraExpressed_names = ','.join(str(e) for e in infraExpressed_names)

    output = {"over":overExpressed_names,"infra":infraExpressed_names,"numOver":overExpressed_number,"numInfra":infraExpressed_number}
    return output


class Results_tutorial(TemplateView):
    template = 'results_tutorial.html'

    def get(self, request):
        
        jobID = 'M353XWU1N37BMNP'

        #Get config file

        configFile = open(settings.MEDIA_ROOT+jobID+'/config.json')
        config = json.load(configFile)

        methods = config['methods']
        jobDir = config['jobDir']
        typeJob = config['typeJob']

        start_date = datetime.date.today()
        end_date = start_date + datetime.timedelta(days=15)
        
        ##All visualizations
        visualization = False

        heatmap = []
        pca = []
        downloads = {}
        distribution = {}
        top10 = {}
        top10FC = {}
        de_groups={}
        features = {}
        rleplots = {}

        i = 0
        j = 0

        Info_Groups = []
        try:
            groupsFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","summary","groups.txt"),'r')
            for line in groupsFile:
                Info_Groups.append(line.strip())

        except:
            pass

        FC_groups = []
        try:
            groupsFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"graphs","summary","combinations.txt"),'r')
            for line in groupsFile:
                FC_groups.append(line.strip())
        except:
            pass


        #InfoGain per group
        infoGain = {}
        downloads['info'] = {}
        for group in Info_Groups:
            info = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".html")
            infoPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".png")
            infoGain[group] = {'HTML':info,'PNG':infoPNG,'name':"Information Gain "}
            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".txt")
            try:
                values = downloads['info'][group]
            except:
                downloads['info'][group] = {}
            downloads['info'][group] = [downloadLink,"infoGain_"+group+".txt"]


        infoGainPairwise = {}
        for group in FC_groups:
            info = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".html")
            infoPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".png")
            infoGainPairwise[group] = {'HTML':info,'PNG':infoPNG,'name':"Information Gain "}

            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","infoGain_"+group+".txt")
            try:
                values = downloads['info'][group]
            except:
                downloads['info'][group] = {}
            downloads['info'][group] = [downloadLink,"infoGain_"+group+".txt"]
        

        ## Features
        no_normalized = os.path.join(settings.MEDIA_ROOT,jobID,"normalized","matrix_NN.txt")
        features = (pd.read_table(no_normalized)['name']).tolist()
        features.sort()

        ##Dif expr
        diffExpr = {}
        consensus = {}
        topDEPerMethod = {}
        downloads['DE']={}
        if config['diffExpr'] == 'True':
            try:
                groupsFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"DE","groups.txt"),'r')
                groups = []
                for line in groupsFile:
                    groups.append(line.strip())
                de_groups = groups
                
                for group in groups:
                    group = group.replace("-","_")
                    groupFile = group
                    group1 = group.split("_")[0]
                    group2 = group.split("_")[1]
                    summaryDE,resultsGroup,consensusM = de_prepare(jobID,group,group1,group2)
                    group = group.replace("_","-")
                    diffExpr[group]= summaryDE,resultsGroup
                    consensus[group] = consensusM
                    topDEPerMethod[group] = de_prepareTop(jobID,group)

                    #Downloads DE
                    protocols = {}
                    protocols = {"edgeR":"edgeR","DESeq2":"deseq2","NOISeq":"noiseq","T-Test":"test"}
                    
                    for protocolName in protocols:
                        protocolFile = protocols[protocolName]
                        downloadLink = os.path.join(settings.MEDIA_URL,jobID,"DE",protocolFile+"_"+groupFile+".txt")
                        try:
                            values = downloads['DE'][protocolName]
                        except:
                            downloads['DE'][protocolName] = {}
                        downloads['DE'][protocolName][group] = [downloadLink,protocolFile+"_"+groupFile+".txt"]
            except:
                pass

        ##Downloads All
        downloadLink = os.path.join(settings.MEDIA_URL,jobID,"results_"+jobID+".zip")
        downloads['all'] = {}
        downloads['all'] = [downloadLink,"results_"+jobID+".zip"]

        downloads['normalized'] = {}
        downloads['top10FC'] = {}
        downloads['top10'] = {}
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

            #Downloads
            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10_"+method+".txt")
            try:
                values = downloads['top10'][METHODS[method]]
            except:
                downloads['top10'][METHODS[method]] = []
            downloads['top10'][METHODS[method]] = [downloadLink,method+".txt"]

            #Top miRNAS FC
            top10FC[method] = {}

            if FC_groups:

                top10FC[method]['name'] = METHODS[method]
                top10FC[method]["comparisons"] = FC_groups
                top10FC[method]["data"] = {}

                for element in FC_groups:
                    topFC = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10FC_"+method+"_"+element+".html")
                    topFCPNG = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10FC_"+method+"_"+element+".png")
                    try:
                        value = top10FC[method]["data"][element]
                    except:
                        top10FC[method]["data"][element] = {}
                    top10FC[method]["data"][element]['HTML'] = topFC
                    top10FC[method]["data"][element]['PNG'] = topFCPNG
                    top10FC[method]["data"][element]['name'] = METHODS[method]
                    if j==0:
                        top10FC[method]["data"][element]['active'] = "block;"
                    else:
                        top10FC[method]["data"][element]['active'] = "none;"

                    #Downloads
                    downloadLink = os.path.join(settings.MEDIA_URL,jobID,"graphs","summary","top10FC_"+method+"_"+element+".txt")
                    try:
                        values = downloads['top10FC'][METHODS[method]]
                    except:
                        downloads['top10FC'][METHODS[method]] = {}
                    downloads['top10FC'][METHODS[method]][element] = [downloadLink,method+"_"+element+".txt"]

                    j = 1
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

                        #RLE plots
            
            if method=="NN":
                rleplotNN = {}
                NN_PNG = os.path.join(settings.MEDIA_URL,jobID,"normalized","RLEplot_NN.png")
                rleplotNN['PNG'] = NN_PNG
                rleplotNN['name'] = METHODS[method]
            else:
                rleplots[method] = {}
                rlePNG = os.path.join(settings.MEDIA_URL,jobID,"normalized","RLEplot_"+method+".png")
                rleplots[method]['PNG'] = rlePNG
                rleplots[method]['name'] = METHODS[method]
                if i==1:
                    rleplots[method]['active'] = "block;"
                else:
                    rleplots[method]['active'] = "none;"

            ##Downloads Normalized
            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"normalized","matrix_"+method+".txt")
            downloads['normalized'][METHODS[method]] = [downloadLink,"matrix_"+method+".txt"]

            i = i + 1 

        summary = {'distribution':distribution,'top10':top10,'top10FC':top10FC,'FC_groups':FC_groups,'info':infoGain,'infoPairwise':infoGainPairwise,'rleplotNN':rleplotNN,'rleplots':rleplots}

        if heatmap or pca:
            visualization=True

        de_software = ["edgeR","DESeq2","NOISeq","TTest"]
        return render(request, self.template, {"jobID":jobID,"typeJob":typeJob,"visualization":visualization,"heatmapPlots":heatmap,
        "pcaPlots":pca,"batchEffect":batchEffect,"downloads":downloads,"summary":summary,"de":diffExpr,"de_groups":de_groups,'date':end_date,
        'consensus':consensus,'de_software':de_software, 'topDEPerMethod':topDEPerMethod,'features':features})
