from django import forms
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, HTML, Field, ButtonHolder,Button
from crispy_forms.bootstrap import Tab,TabHolder,Div
from django.core.validators import MaxValueValidator, MinValueValidator 

from query.utils import create_card,render_modal
from normSeq.settings import *


class Query(forms.Form):
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


    matrix = forms.FileField(label='Upload input matrix', required=False)
    url = forms.URLField(label='URL/link', required=False)
    jobID = forms.CharField(widget = forms.HiddenInput())
    typeJob = forms.CharField(widget = forms.HiddenInput(), initial = "miRNA")
    methods = forms.MultipleChoiceField(choices=METHODS,label="Choose one or more normalization methods:",required=True)
    mirgenedb = forms.ChoiceField(choices=METHODS,label="Choose short name from MiRGeneDB:",required=False)
    annotation = forms.FileField(label="Annotation:", required=False)
    annotationURL = forms.URLField(label='URL/link to annotation file:', required=False)
    minrc = forms.IntegerField(label="Minimum RC",validators=[MinValueValidator(0)],initial=0,required=False)
    
    batchEffect = forms.TypedChoiceField(coerce=lambda x: x =='True', 
                                   choices=((False, 'No'), (True, 'Yes')),label="Apply batch effect correction:")
    batchFile = forms.FileField(label="Batch Effect Annotation File:", required=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.layout = Layout(
            create_card(
                TabHolder(
                        Tab(' Upload ',
                            HTML("""<br>"""),
                            Field('matrix', css_class='form-control'),
                            css_id='upload'
                        ),
                        Tab(' URL ',
                            HTML("""<br>"""),
                            Field('url', css_class='form-control',placeholder="e.g.:https://url.com/matrix.txt"),
                            css_class="link",
                            css_id='url'),
                        css_id='tab_upload'
                ),
                title='Input Matrix',
                id="input",
                show=True,
                modal="Choose_Input"
            ),
            HTML("""<br>"""),
            create_card(
                TabHolder(
                    Tab('Upload Annotation File',
                    HTML("""<br>"""),
                    Field('annotation',css_class='form-control'),
                    css_id='annotation_upload'
                    ),
                    Tab('URL of Annotation File',
                    HTML("""<br>"""),
                    Field('annotationURL',css_class='form-control'),
                    css_id = 'annotation_url'
                    ),
                    css_id="tab_annotation"
                ),
                HTML('<a href="'+STATIC_URL+'testFiles/template.txt" download="template.txt"><button type="button" class="btn btn-info btn-sm float-right">Download sample template</button></a>'),
                title="Annotation",
                id="annotation",
                show=True,
                modal="Choose_Annotation"
            ),
            HTML("""<br>"""),
            create_card(
                Field('minrc',css_class='form-control'),
                Field('batchEffect',css_class='form-control',css_id='batchE'),
                Div(
                    Field('batchFile', css_class='form-control'),
                    HTML('<a href="'+STATIC_URL+'testFiles/template_batch.txt" download="template.txt"><button type="button" class="btn btn-info btn-sm float-right mb-3">Download sample template</button></a>'),
                    css_id='divFileBatch'),
                title="Parameters",
                id="parameters",
                show=False,
                modal="chooseParameters"
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
                Submit('submit', 'SUBMIT', css_class='btn btn-success btn-xl', css_id="submit", onclick='check_form();')
            )
        )
    def clean(self):
        ALLOWED_FORMATS = ["txt","xlsx","tsv","csv"]
        cleaned_data = super(Query, self).clean()

        if not cleaned_data.get('matrix') and not cleaned_data.get('url'):
             self.add_error('matrix', 'One of these two fields is required (file or url)')
             self.add_error('url', 'One of these two fields is required (file or url)')
        if cleaned_data.get('matrix') and cleaned_data.get('url'):
            self.add_error('matrix', 'Choose either file or URL')
            self.add_error('url', 'Choose either file or URL')

        if cleaned_data.get('matrix'):
            fileType = cleaned_data.get('matrix').name.split(".")[1]
            if not fileType in ALLOWED_FORMATS:
                self.add_error('matrix', fileType+' is not an allowed matrix format. Please choose another file.')
        

        return cleaned_data
