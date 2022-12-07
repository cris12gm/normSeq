from django import forms
from django.core.validators import RegexValidator
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, Div,HTML,Field, ButtonHolder
from crispy_forms.bootstrap import (
    Tab,
    TabHolder,
    Accordion,
    AccordionGroup
)

from query.utils import create_card,render_modal




class QueryMirna(forms.Form):
    METHODS = (
        ('NN', ' No normalization, just visualization'),
        ('CPM', ' Counts per Million'),
        ('TC', ' Total Count'),
        ('UQ', ' Upper Quartile'),
        ('Med',' Median'),
        ('TMM',' TMM'),
        ('QN',' Quantile'),
#        ('RUV',' Remove Unwanted Variation'),
        ('RLE','Relative Log Expression')
    )


    matrix = forms.FileField(required=False)
    url = forms.CharField(required=False)
    jobID = forms.CharField(widget = forms.HiddenInput())
    typeJob = forms.CharField(widget = forms.HiddenInput(), initial = "miRNA")
    methods = forms.MultipleChoiceField(choices=METHODS,label="Choose one or more normalization methods:",required=True)
    mirgenedb = forms.ChoiceField(choices=METHODS,label="Choose short name from MiRGeneDB:",required=False)
    annotation = forms.FileField(label="Annotation:", required=False)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.layout = Layout(
            create_card(
                TabHolder(
                        Tab(' Upload ',
                            HTML("""<br>"""),
                            Field('matrix', css_class='form-control'),
                        ),
                        Tab(' URL ',
                            HTML("""<br>"""),
                            Field('url', css_class='form-control',placeholder="e.g.:https://url.com/matrix.txt"),
                            css_class="link")
                ),
                title='Input Matrix',
                id="input",
                show=True,
                modal="Choose_Input"
            ),
            HTML("""<br>"""),
            create_card(
                Field('annotation',css_class='form-control'),
                HTML("""<br>"""),
                Field('mirgenedb',css_class='form-control'),
                title="Annotation",
                id="annotation",
                show=False,
                modal="Choose_Annotation"
            ),
            HTML("""<br>"""),
            create_card(
                Field('methods'),
                title= "Normalization Methods",
                id="normalization",
                show=True,
                modal="Choose_Method"
            ),
            HTML("""<br>"""),
            'jobID',
            'typeJob',
            ButtonHolder(
                Submit('submit', 'SUBMIT', css_class='btn btn-success btn-xl', onclick='check_form();')
            )
        )