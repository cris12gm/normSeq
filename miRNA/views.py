from django.conf import settings
import random,string
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic.edit import CreateView
from django.views.generic import FormView, DetailView, TemplateView
from django.http import JsonResponse
from django.urls import reverse_lazy
from django.core.files.storage import FileSystemStorage
from .forms import Query
from query.models import Job
from .utils import createGridPlot


import os
import pandas as pd
import json
import math

import time

class Errors(Enum):
    NO_ERROR = 0
    NOT_VALID = 1
    NOT_ASSOCIATED = 2

class miRNAResults(TemplateView):
    template = 'mirna.html'

    def get(self, request):

        if not request.GET['jobID']:
            return redirect(settings.SUB_SITE+"/query/")

        jobID = request.GET['jobID']

        #Get config file

        configFile = open(settings.MEDIA_ROOT+jobID+'/config.json')
        config = json.load(configFile)

        methods = config['methods']
        jobDir = config['jobDir']

        ##All visualizations
        visualization = False

        heatmap = []

        for method in methods:
            fileHeatmap = os.path.join(jobDir,"graphs","heatmap_"+method+".html")
            with open(fileHeatmap,'r') as file:
                heatmapHTML = file.read().rstrip()
            pngHeatmap = os.path.join(settings.MEDIA_URL,jobID,"graphs","heatmap_"+method+".png")
            id_modal = "heatmap_"+method
            heatmap.append([pngHeatmap,heatmapHTML,id_modal,method])
        
        visualization=True

        return render(request, self.template, {"jobID":jobID,"visualization":visualization,"heatmapPlots":heatmap})

