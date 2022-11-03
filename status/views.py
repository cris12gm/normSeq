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
from query.models import Job

class Errors(Enum):
    NO_ERROR = 0
    NOT_VALID = 1
    NOT_ASSOCIATED = 2

class statusJob(TemplateView):
    template = 'status.html'

    jobID = ""
    def get(self, request):  

        return redirect(settings.SUB_SITE+"/query/")

    def post(self, request):
        jobID = request.POST['jobID']

        jobDB = Job.objects.get(JobID=jobID)
        status = jobDB.getStatus()

        return render(request, self.template, {"jobID":jobID,"status":status})
    
    