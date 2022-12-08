from django import forms
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, HTML, Field, ButtonHolder,Button
from crispy_forms.bootstrap import Tab,TabHolder

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




    matrix = forms.FileField(label='Upload input matrix', required=False)
    url = forms.URLField(label='URL/link', required=False)
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
                HTML('<button type="button" class="btn btn-info btn-sm float-right">Download sample template</button>'),
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
    def clean(self):
        ALLOWED_FORMATS = ["txt","xlsx","tsv","csv"]
        cleaned_data = super(QueryMirna, self).clean()


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
