from django.conf import settings
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic import TemplateView
from query.models import Job

import os
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

        try:
            jobDB = Job.objects.get(JobID=jobID)
        except:
            return redirect(settings.SUB_SITE+"/query/")
        status = jobDB.getStatus()
        typeJob = jobDB.getType()

        #Check if the job has finished
        pid = jobDB.getPid()
        resultsFile = os.path.join(settings.BASE_DIR,settings.MEDIA_ROOT+jobID,"results.txt")
        statusFile = os.path.join(settings.BASE_DIR,settings.MEDIA_ROOT+jobID,"status.txt")
        errorFile = os.path.join(settings.BASE_DIR,settings.MEDIA_ROOT+jobID,"Error.txt")
        
        messages = {}
        try:
            messages = open(statusFile,'r').readlines()[0]
        except:
            pass
        try:
            p = psutil.Process(pid)
            statusPid = (p.status() if hasattr(p.status, '__call__'
                                        ) else p.status)
            if statusPid=="zombie" and os.path.exists(resultsFile):
                jobDB.alterStatus("Finished")

        except:
            if os.path.exists(resultsFile):
                jobDB.alterStatus("Finished")
            else:
                jobDB.alterStatus("Error")
        if status=="Created":
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob})
        elif status=="Finished":
            return redirect(settings.SUB_SITE+"/results/?jobID="+jobID)
        elif status=="Running":
            configFile = jobDB.getConfig()
            config = json.load(open(configFile,'r'))
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob,"messages":messages,"config":config})
        else:
            try:
                messages = open(errorFile,'r').readlines()[0]
            except:
                messages = "There is an error with your job, please report it <a target='_blank' href='https://github.com/cris12gm/normSeq/issues'>here</a> with the following ID <b>"+jobID+"</b>"
            return render(request, self.template, {"jobID":jobID,"status":status,"error":messages,"typeJob":typeJob})
            

    def post(self, request):
        jobID = request.POST['jobID']

        #Check if jobID exists
        try:
            jobDB = Job.objects.get(JobID=jobID)
            return redirect(settings.SUB_SITE+"/status?jobID="+jobID)
        except:
            return redirect(settings.SUB_SITE+"/query")

        
    

    