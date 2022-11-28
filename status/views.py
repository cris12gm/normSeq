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

import os
import subprocess
import psutil
import json

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

        #Check if the job has finished
        pid = jobDB.getPid()
        resultsFile = os.path.join(settings.BASE_DIR,settings.MEDIA_ROOT+jobID,"results.txt")

        if not psutil.pid_exists(pid) and os.path.exists(resultsFile):
            jobDB.alterStatus("Finished")
        elif not psutil.pid_exists(pid) and not os.path.exists(resultsFile):
            jobDB.alterStatus("Error")



        if status=="Created":
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob})
        elif status=="Finished":
            if typeJob=="miRNA":
                return redirect(settings.SUB_SITE+"/mirna/?jobID="+jobID)
        elif status=="Running":
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob})
        else:
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob}) ##Poner error aqui
            

    def post(self, request):
        jobID = request.POST['jobID']
        typeJobForm = request.POST['typeJob']
        methodsForm = request.POST.getlist('methods')

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

            typeJob = jobDB.getType()

            #Parameters contain the info for config.json
            parameters = {}
            parameters["jobID"]=jobID
            parameters["typeJob"]=typeJob
            parameters["methods"]=methodsForm

            #Launch
            if typeJob =="miRNA":
                
                #Save the rest of parameters and launch

                with open(settings.MEDIA_ROOT+jobID+'/config.json', 'w') as fp:
                    json.dump(parameters, fp)

                script = os.path.join(settings.BASE_DIR,"bin","miRNA_bench.py")
                
                jobDir = os.path.join(settings.BASE_DIR,settings.MEDIA_ROOT+jobID)
                scriptCm = ["python3",script,jobDir]
                #scriptCm = ["python3",script,jobDir]
                #response = subprocess.Popen(' '.join(scriptCm),shell=True).pid
                response = subprocess.Popen(' '.join(scriptCm)).pid
                if psutil.pid_exists(response):
                    jobDB.alterPid(response)
                    jobDB.alterStatus("Running")
                else:
                    jobDB.alterStatus("Error")
        
        return redirect(settings.SUB_SITE+"/status?jobID="+jobID)
    

    