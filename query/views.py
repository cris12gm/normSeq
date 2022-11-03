from django.conf import settings
import random,string
from enum import Enum
from functools import reduce
from django.shortcuts import render, redirect
from django.views.generic.edit import CreateView
from django.views.generic import FormView, DetailView, TemplateView
from django.http import JsonResponse
from django.urls import reverse_lazy
from .forms import QueryMirna
from .models import Job


def generate_uniq_id(size=15, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def generate_id():
    is_new = True
    while is_new:
        pipeline_id = generate_uniq_id()
        #if not JobStatus.objects.filter(pipeline_key=pipeline_id):
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
            Job.objects.create(JobID= jobID, status="Created")
        formMirna = QueryMirna()

        return render(request, self.template, {"jobID":jobID,"formMirna":formMirna,
        })
    