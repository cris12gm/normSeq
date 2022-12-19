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
import math

import time


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

        i = 0
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
             
            ##Downloads
            downloadLink = os.path.join(settings.MEDIA_URL,jobID,"normalized","matrix_"+method+".txt")

            downloads[METHODS[method]] = [downloadLink,"matrix_"+method+".txt"]
            i = i + 1 

        summary = {'distribution':distribution,'top10':top10}

        if heatmap or pca:
            visualization=True

        return render(request, self.template, {"jobID":jobID,"typeJob":typeJob,"visualization":visualization,"heatmapPlots":heatmap,
        "pcaPlots":pca,"downloads":downloads,"summary":summary})


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