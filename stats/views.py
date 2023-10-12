from django.conf import settings
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic import TemplateView
from query.models import Job


class Errors(Enum):
    NO_ERROR = 0
    NOT_VALID = 1
    NOT_ASSOCIATED = 2

class statisticsDB(TemplateView):
    template = 'statistics.html'

    def get(self, request):  
        
        success,error = Job.stats()
        
        return render(request, self.template, {'success':success,'error':error})
            

        
    

    