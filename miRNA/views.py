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

class Errors(Enum):
    NO_ERROR = 0
    NOT_VALID = 1
    NOT_ASSOCIATED = 2

class miRNAResults(TemplateView):
    template = 'mirna.html'
    jobID = ""
    def get(self, request):
        print(request.GET)
        try:
            jobID = request.GET['jobID']
        except:
            jobID = ""
        return render(request, self.template, {"jobID":jobID})