from django import forms
from django.core.validators import RegexValidator

my_validator = RegexValidator(r"rs[0-9]{1,15}", "The snpID should be in this format: rs727563")
class Query(forms.Form):
    SNPid = forms.CharField(label='SNP ID',max_length=15, validators=[my_validator], widget=forms.TextInput())
    #attrs={'placeholder':'rs727563',} -> inside TextInput
