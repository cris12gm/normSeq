from django.db import models
from django import forms
from django.utils import timezone


class Job(models.Model):
    JobID = models.CharField(max_length=15)
    status = models.CharField(max_length=200)
    typeJob = models.CharField(max_length=25)
    created_date = models.DateTimeField(
            default=timezone.now)


    def __str__(self):
        return self.JobID

    def getStatus(self):
        return self.status

    def getType(self):
        return self.typeJob
    
    def alterType(self,typeJobForm):
        self.typeJob = typeJobForm
        self.save()
