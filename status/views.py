from django.conf import settings
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic import TemplateView
from query.models import Job

import os
import psutil

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
            if typeJob=="miRNA":
                return redirect(settings.SUB_SITE+"/mirna/?jobID="+jobID)
        elif status=="Running":
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob})
        else:
            return render(request, self.template, {"jobID":jobID,"status":status,"typeJob":typeJob}) ##Poner error aqui
            

    def post(self, request):
        jobID = request.POST['jobID']

        return redirect(settings.SUB_SITE+"/status?jobID="+jobID)
    

    