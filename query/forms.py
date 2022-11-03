from django import forms
from django.core.validators import RegexValidator



class QueryMirna(forms.Form):

    matrix = forms.FileField()
    jobID = forms.CharField(widget = forms.HiddenInput())