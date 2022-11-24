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


import os
import pandas as pd

import time

class Errors(Enum):
    NO_ERROR = 0
    NOT_VALID = 1
    NOT_ASSOCIATED = 2

class miRNAResults(TemplateView):
    template = 'mirna.html'

    def get(self, request):

            ##All visualizations
        visualization = True
        heatmapPlot = False

        if not request.GET['jobID']:
            return redirect(settings.SUB_SITE+"/query/")

        jobID = request.GET['jobID']
        
        rpmTable = pd.read_table(settings.MEDIA_ROOT+jobID+"/matrix_RPM.txt",sep="\t",index_col="name")
        
        #Visualizations
            #Heatmap

        
        
        fileHeatmap = settings.MEDIA_ROOT+jobID+"/heatmap.html"

        if os.path.exists(fileHeatmap):
            with open(fileHeatmap,'r') as file:
                heatmapPlot = file.read().rstrip()
                visualization = "True"

        return render(request, self.template, {"jobID":jobID,"rpmTable":rpmTable,"visualization":visualization,"heatmapPlot":heatmapPlot})