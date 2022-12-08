from django.conf import settings
import random,string
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic import  TemplateView
from .forms import QueryMirna
from .models import Job
from django.core.files.storage import FileSystemStorage

import os,sys
import subprocess
import psutil
import json

def generate_uniq_id(size=15, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def generate_id():
    is_new = True
    while is_new:
        pipeline_id = generate_uniq_id()
        if not Job.objects.filter(JobID=pipeline_id):
            return pipeline_id

class Errors(Enum):
    NO_ERROR = 0
    NOT_VALID = 1
    NOT_ASSOCIATED = 2

class Query(TemplateView):
    template = 'query.html'

    def get(self, request):  

        jobID = generate_id()
        if jobID:
            Job.objects.create(JobID= jobID, status="Created", typeJob="Created")
        formMirna = QueryMirna()

        return render(request, self.template, {"jobID":jobID,"formMirna":formMirna,
        
        })

    def post(self,request):

        form = QueryMirna(request.POST,request.FILES)

        if form.is_valid():
            #Validate post and redirect
            
            jobID = request.POST['jobID']
            typeJobForm = request.POST['typeJob']
            methodsForm = request.POST.getlist('methods')

            #Get object from models
        
            jobDB = Job.objects.get(JobID=jobID)
            typeJob = jobDB.getType()

            if typeJob=="Created":
            #Alter type to match form selected
                jobDB.alterType(typeJobForm)
                typeJob = jobDB.getType()
            
            #Save file
            try:
                matrixFile = request.FILES['matrix']
                extension = matrixFile.name.split(".")[1]  
                fs = FileSystemStorage()
                filename = fs.save(jobID+"/"+"matrix."+extension, matrixFile)
                uploaded_file_url = fs.url(filename)
            except:
                try:
                    url = request.FILES['url']
                except:
                    pass

            typeJob = jobDB.getType()


            #Save annotation

            try:
                annotationFile = request.FILES['annotation']
                annotation = annotationFile.name
                fs = FileSystemStorage()
                filename = fs.save(jobID+"/"+annotation, annotationFile)
            except:
                annotation=False

            # Parameters contain the info for config.json
            parameters = {}
            parameters["inputExtension"] = extension
            parameters["jobID"] = jobID
            parameters["typeJob"] = typeJob
            parameters["methods"] = methodsForm
            parameters["annotation"] = annotation
            parameters["jobDir"]= os.path.join(settings.MEDIA_ROOT+jobID)

            # Launch
            if typeJob =="miRNA":
                
                #Save the rest of parameters and launch

                with open(settings.MEDIA_ROOT+jobID+'/config.json', 'w') as fp:
                    json.dump(parameters, fp)

                script = os.path.join(settings.BASE_DIR,"bin","miRNA_bench.py")
                
                
                configFile = settings.MEDIA_ROOT+jobID+'/config.json'

                scriptCm = ["python",script,configFile]

                response = subprocess.Popen(' '.join(scriptCm),shell=True).pid
                if psutil.pid_exists(response):
                    jobDB.alterPid(response)
                    jobDB.alterStatus("Running")
                else:
                    jobDB.alterStatus("Error")

            return redirect(settings.SUB_SITE+"/status?jobID="+jobID)
        else:
            jobID = request.POST['jobID']
            form = QueryMirna(request.POST,request.FILES)
        return render(request, self.template, {"formMirna":form,"jobID":jobID})