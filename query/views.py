from django.conf import settings
import random,string
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic import  TemplateView
from .forms import QueryMirna,QueryMrna,QueryOther,QueryTrna
from .models import Job
from django.core.files.storage import FileSystemStorage

import os,sys
import pandas as pd
import subprocess
import psutil
import json
import requests

ALLOWED_FORMATS = ["txt","xlsx","tsv","csv"]

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
        formTrna = QueryTrna()
        formMrna = QueryMrna()
        formOther = QueryOther()

        return render(request, self.template, {"jobID":jobID,"formMirna":formMirna,"formMrna":formMrna,
        "formTrna":formTrna,"formOther":formOther})

    def post(self,request):

        rnaClass = request.POST['typeJob']

        if rnaClass=="miRNA":
            form = QueryMirna(request.POST,request.FILES)
            if form.is_valid():
                #Validate post and redirect
            
                jobID = request.POST['jobID']
                typeJobForm = request.POST['typeJob']
                methodsForm = request.POST.getlist('methods_mirna')

                # Get object from models
                jobDB = Job.objects.get(JobID=jobID)
                typeJob = jobDB.getType()

                if typeJob=="Created":
            
                #Alter type to match form selected
                    jobDB.alterType(typeJobForm)
                    typeJob = jobDB.getType()
            
                try:
                    matrixFile = request.FILES['matrix_mirna']

                    extension = matrixFile.name.split(".")[1]

                    fs = FileSystemStorage()
                    filename = fs.save(jobID+"/"+"matrix."+extension, matrixFile)

                    if extension == "csv":
                        df = pd.read_csv(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.csv"),sep=";")
                        df2=df.dropna(how='all')
                        df2.to_csv(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.txt"),sep="\t",index=None)
                    elif extension == "tsv":
                        filename = fs.save(jobID+"/"+"matrix.txt", matrixFile)
                    elif extension == "xlsx":
                        df = pd.read_excel(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.xlsx"))
                        df2=df.dropna(how='all')
                        df2.to_csv(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.txt"),sep="\t",index=None)
                except:
                    try:
                        url = request.POST['url_mirna']
                        filename = url.split("/")[-1]
                        extension = filename.split(".")[1]
                        response = requests.get(url)
                        contenido = response.content.decode("utf-8")
                        if not os.path.exists(os.path.join(settings.MEDIA_ROOT,jobID)):
                            os.mkdir(os.path.join(settings.MEDIA_ROOT,jobID))
                        matrixFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"matrix."+extension),'a')
                        matrixFile.write(contenido)
                        matrixFile.close()

                        if extension == "csv":
                            df = pd.read_csv(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.csv"),sep=";")
                            df2=df.dropna(how='all')
                            df2.to_csv(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.txt"),sep="\t",index=None)
                        elif extension == "tsv":
                            filename = fs.save(jobID+"/"+"matrix.txt", matrixFile)
                        elif extension == "xlsx":
                            df = pd.read_excel(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.xlsx"))
                            df2=df.dropna(how='all')
                            df2.to_csv(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.txt"),sep="\t",index=None)

                    except:
                        return render(request, "error.html", {})
                        #Error

                #Save annotation

                try:
                    try:
                        annotationFile = request.FILES['annotation_mirna']
                        annotation = annotationFile.name
                        fs = FileSystemStorage()
                        filename = fs.save(jobID+"/"+"annotation.txt", annotationFile)
                        annotation=True
                    except:
                        annotationURL = request.POST['annotationURL_mirna']
                        filename = annotationURL.split("/")[-1]
                        response = requests.get(annotationURL)
                        contenido = response.content.decode("utf-8")
                        matrixFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"annotation.txt"),'a')
                        matrixFile.write(contenido)
                        matrixFile.close()
                        annotation = True

                except:
                    #Read annotation from matrix
                    samples = open(os.path.join(settings.MEDIA_ROOT,jobID,"matrix.txt")).readline()
                    samples = samples.strip().split("\t")
                    samples.pop(0)
                    annotationFile = open(os.path.join(settings.MEDIA_ROOT,jobID,"annotation.txt"),'a')
                    annotationFile.write("sample\tgroup\n")
                    for sample in samples:
                        annotationFile.write(sample+"\t"+sample+"\n")
                    annotationFile.close()

                    annotation=False

                # Parameters contain the info for config.json
                parameters = {}
                parameters["inputExtension"] = extension
                parameters["jobID"] = jobID
                parameters["typeJob"] = typeJob
                parameters["methods"] = methodsForm
                parameters["annotation"] = annotation
                parameters["jobDir"]= os.path.join(settings.MEDIA_ROOT+jobID)

                # Save parameters and Launch
           
           
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
                formMirna = QueryMirna(request.POST,request.FILES)
                formTrna = QueryTrna()
                formMrna = QueryMrna()
                formOther = QueryOther()

                return render(request, self.template, {"formMirna":formMirna,"jobID":jobID,
                "formMrna":formMrna,"formTrna":formTrna,"formOther":formOther})
        else:
            jobID = request.POST['jobID']
            formMirna = QueryMirna(request.POST,request.FILES)
            formTrna = QueryTrna()
            formMrna = QueryMrna()
            formOther = QueryOther()

            return render(request, self.template, {"formMirna":formMirna,"jobID":jobID,
                "formMrna":formMrna,"formTrna":formTrna,"formOther":formOther})