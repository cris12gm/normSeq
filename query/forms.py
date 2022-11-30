from django import forms
from django.core.validators import RegexValidator
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, Div,HTML,Field
from crispy_forms.bootstrap import TabHolder, Tab
from query.utils import create_card,render_modal

class QueryMirna(forms.Form):
    METHODS = (
        ('CPM', ' Counts per Million'),
        ('TC', ' Total Count'),
        ('UQ', ' Upper Quartile'),
        ('Med',' Median'),
        ('DESEQ',' DEseq'),
        ('TMM',' TMM (edgeR)'),
        ('QN',' Quantile'),
        ('RUV',' Remove Unwanted Variation'),
        ('NN', ' No normalization')
    )
    
    matrix = forms.FileField()
    url = forms.CharField()
    jobID = forms.CharField(widget = forms.HiddenInput())
    typeJob = forms.CharField(widget = forms.HiddenInput(), initial = "miRNA")
    methods = forms.MultipleChoiceField(choices=METHODS,label="Choose one or more normalization methods:",required=True)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.layout = Layout(
            create_card(
                TabHolder(
                    Tab(' Upload ',
                        HTML("""<br>"""),
                        Field('matrix', css_class='form-control'),
                    ),
                    #Tab(' URL/link ',
                    #    HTML("""<br>"""),
                    #    Field('url', css_class='form-control'),
                    #)
                ),
                title="Input"+render_modal("Choose_Input")
            ),
            HTML("""<br>"""),
            Div(
                create_card(
                    Div(Field('methods')),title="Methods"+render_modal("Choose_Method"))
            ),
            'jobID',
            'typeJob',
            HTML("""<br>"""),
            Submit('submit', 'Submit', css_class='btn btn-success btn-sm')
        )
