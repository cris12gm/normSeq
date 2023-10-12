from django.db import models
from django import forms
from django.utils import timezone


class Job(models.Model):
    JobID = models.CharField(max_length=15)
    status = models.CharField(max_length=200)
    typeJob = models.CharField(max_length=25)
    configFile = models.CharField(max_length=1000)
    pid = models.IntegerField(default=0)
    created_date = models.DateTimeField(
            default=timezone.now)


    def __str__(self):
        return self.JobID

    def getStatus(self):
        return self.status

    def getType(self):
        return self.typeJob
    
    def getPid(self):
        return self.pid
    
    def getConfig(self):
        return self.configFile
    
    def getDate(self):
        return self.created_date

    def alterType(self,typeJobForm):
        self.typeJob = typeJobForm
        self.save()

    def alterStatus(self,newStatus):
        self.status = newStatus
        self.save()

    def alterPid(self,newPID):
        self.pid = newPID
        self.save()

    def stats():
        success = {}
        success['miRNA'] = Job.objects.filter(status='Finished',typeJob='miRNA').count()
        success['mRNA'] =  Job.objects.filter(status='Finished',typeJob='mRNA').count()
        success['tRNA'] =  Job.objects.filter(status='Finished',typeJob='tRNA').count()
        success['other'] = Job.objects.filter(status='Finished',typeJob='Feature').count()
        success['total'] = Job.objects.filter(status='Finished').count()

        error = Job.objects.filter(status='Error').count()

        return success,error