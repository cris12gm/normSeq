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

    def get(self, request):  

        try:
            jobID = request.GET['jobID']
        except:
            return redirect(settings.SUB_SITE+"/query/")

        jobDB = Job.objects.get(JobID=jobID)
        status = jobDB.getStatus()
        typeJob = jobDB.getType()

        if status=="Created":
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob})
        elif status=="Finished":
            if typeJob=="miRNA":
 #               return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob})
                return redirect(settings.SUB_SITE+"/mirna/?jobID="+jobID)
            

    def post(self, request):
        jobID = request.POST['jobID']
        typeJobForm = request.POST['typeJob']

        #Get object from models
        jobDB = Job.objects.get(JobID=jobID)
        typeJob = jobDB.getType()

        # Alter type to match form selected
        if typeJob=="Created":
            jobDB.alterType(typeJobForm)

            #Save file
            matrixFile = request.FILES['matrix']    
            fs = FileSystemStorage()
            filename = fs.save(jobID+"/"+"matrix.txt", matrixFile)
            uploaded_file_url = fs.url(filename)
        return redirect(settings.SUB_SITE+"/status?jobID="+jobID)
    

    