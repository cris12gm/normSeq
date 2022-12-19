from django import forms
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, HTML, Field, ButtonHolder,Button
from crispy_forms.bootstrap import Tab,TabHolder,Div
from django.core.validators import MaxValueValidator, MinValueValidator 

from query.utils import create_card,render_modal
from normSeq.settings import *


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


    matrix_mirna = forms.FileField(label='Upload input matrix', required=False)
    url_mirna = forms.URLField(label='URL/link', required=False)
    jobID = forms.CharField(widget = forms.HiddenInput())
    typeJob = forms.CharField(widget = forms.HiddenInput(), initial = "miRNA")
    methods_mirna = forms.MultipleChoiceField(choices=METHODS,label="Choose one or more normalization methods:",required=True)
    mirgenedb_mirna = forms.ChoiceField(choices=METHODS,label="Choose short name from MiRGeneDB:",required=False)
    annotation_mirna = forms.FileField(label="Annotation:", required=False)
    annotationURL_mirna = forms.URLField(label='URL/link to annotation file:', required=False)
    minrc_mirna = forms.IntegerField(label="Minimum RC",validators=[MinValueValidator(0)],initial=0,required=False)
    
    batchEffect_mirna = forms.TypedChoiceField(coerce=lambda x: x =='True', 
                                   choices=((False, 'No'), (True, 'Yes')),label="Apply batch effect correction:")
    batchFile_mirna = forms.FileField(label="Batch Effect Annotation File:", required=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.layout = Layout(
            create_card(
                TabHolder(
                        Tab(' Upload ',
                            HTML("""<br>"""),
                            Field('matrix_mirna', css_class='form-control'),
                            css_id='upload_mirna'
                        ),
                        Tab(' URL ',
                            HTML("""<br>"""),
                            Field('url_mirna', css_class='form-control',placeholder="e.g.:https://url.com/matrix.txt"),
                            css_class="link",
                            css_id='url_mirna'),
                        css_id='tab_upload_mirna'
                ),
                title='Input Matrix',
                id="input_mirna",
                show=True,
                modal="Choose_Input"
            ),
            HTML("""<br>"""),
            create_card(
                TabHolder(
                    Tab('Upload Annotation File',
                    HTML("""<br>"""),
                    Field('annotation_mirna',css_class='form-control'),
                    css_id='annotation_upload_mirna'
                    ),
                    Tab('URL of Annotation File',
                    HTML("""<br>"""),
                    Field('annotationURL_mirna',css_class='form-control'),
                    css_id = 'annotation_url_mirna'
                    ),
                    css_id="tab_annotation_mirna"
                ),
                HTML('<a href="'+STATIC_URL+'testFiles/template.txt" download="template.txt"><button type="button" class="btn btn-info btn-sm float-right">Download sample template</button></a>'),
                title="Annotation",
                id="annotation_mirna",
                show=True,
                modal="Choose_Annotation"
            ),
            HTML("""<br>"""),
            create_card(
                Field('minrc_mirna',css_class='form-control'),
                Field('batchEffect_mirna',css_class='form-control',css_id='batchE'),
                Div(
                    Field('batchFile_mirna', css_class='form-control'),
                    HTML('<a href="'+STATIC_URL+'testFiles/template.txt" download="template.txt"><button type="button" class="btn btn-info btn-sm float-right mb-3">Download sample template</button></a>'),
                    css_id='divFileBatch'),
                title="Parameters",
                id="parameters_mirna",
                show=False,
                modal="chooseParameters"
            ),
            HTML("""<br>"""),
            create_card(
                Field('methods_mirna'),
                title= "Normalization Methods",
                id="normalization_mirna",
                show=True,
                modal="Choose_Method"
            ),
            HTML("""<br>"""),
            'jobID',
            'typeJob',
            ButtonHolder(
                Submit('submit_mirna', 'SUBMIT', css_class='btn btn-success btn-xl', css_id="submit_mirna", onclick='check_form();')
            )
        )
    def clean(self):
        ALLOWED_FORMATS = ["txt","xlsx","tsv","csv"]
        cleaned_data = super(QueryMirna, self).clean()

        if not cleaned_data.get('matrix_mirna') and not cleaned_data.get('url_mirna'):
             self.add_error('matrix_mirna', 'One of these two fields is required (file or url)')
             self.add_error('url_mirna', 'One of these two fields is required (file or url)')
        if cleaned_data.get('matrix_mirna') and cleaned_data.get('url_mirna'):
            self.add_error('matrix_mirna', 'Choose either file or URL')
            self.add_error('url_mirna', 'Choose either file or URL')

        if cleaned_data.get('matrix_mirna'):
            fileType = cleaned_data.get('matrix_mirna').name.split(".")[1]
            if not fileType in ALLOWED_FORMATS:
                self.add_error('matrix_mirna', fileType+' is not an allowed matrix format. Please choose another file.')
        

        return cleaned_data


class QueryMrna(forms.Form):
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


    matrix_mrna = forms.FileField(label='Upload input matrix', required=False)
    url_mrna = forms.URLField(label='URL/link', required=False)
    jobID = forms.CharField(widget = forms.HiddenInput())
    typeJob = forms.CharField(widget = forms.HiddenInput(), initial = "miRNA")
    methods_mrna = forms.MultipleChoiceField(choices=METHODS,label="Choose one or more normalization methods:",required=True)
    mirgenedb_mrna = forms.ChoiceField(choices=METHODS,label="Choose short name from MiRGeneDB:",required=False)
    annotation_mrna = forms.FileField(label="Annotation:", required=False)
    annotationURL_mrna = forms.URLField(label='URL/link to annotation file:', required=False)
    minrc_mrna = forms.IntegerField(label="Minimum RC",validators=[MinValueValidator(0)],initial=0,required=False)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.layout = Layout(
            create_card(
                TabHolder(
                        Tab(' Upload ',
                            HTML("""<br>"""),
                            Field('matrix_mrna', css_class='form-control'),
                            css_id='upload_mrna'
                        ),
                        Tab(' URL ',
                            HTML("""<br>"""),
                            Field('url_mrna', css_class='form-control',placeholder="e.g.:https://url.com/matrix.txt"),
                            css_class="link",
                            css_id='url_mrna'),
                        css_id='tab_upload_mrna'
                ),
                title='Input Matrix',
                id="input_mrna",
                show=True,
                modal="Choose_Input"
            ),
            HTML("""<br>"""),
            create_card(
                TabHolder(
                    Tab('Upload Annotation File',
                    HTML("""<br>"""),
                    Field('annotation_mrna',css_class='form-control'),
                    css_id='annotation_upload_mrna'
                    ),
                    Tab('URL of Annotation File',
                    HTML("""<br>"""),
                    Field('annotationURL_mrna',css_class='form-control'),
                    css_id = 'annotation_url_mrna'
                    ),
                    css_id="tab_annotation_mrna"
                ),
                HTML('<a href="'+STATIC_URL+'testFiles/template.txt" download="template.txt"><button type="button" class="btn btn-info btn-sm float-right">Download sample template</button></a>'),
                title="Annotation",
                id="annotation_mrna",
                show=True,
                modal="Choose_Annotation"
            ),
            HTML("""<br>"""),
            create_card(
                Field('minrc_mrna',css_class='form-control'),
                title="Parameters",
                id="parameters_mrna",
                show=False,
                modal="Parameters"
            ),
            HTML("""<br>"""),
            create_card(
                Field('methods_mrna'),
                title= "Normalization Methods",
                id="normalization_mrna",
                show=True,
                modal="Choose_Method"
            ),
            HTML("""<br>"""),
            'jobID',
            'typeJob',
            ButtonHolder(
                Submit('submit_mrna', 'SUBMIT', css_class='btn btn-success btn-xl', onclick='check_form();')
            )
        )
    def clean(self):
        ALLOWED_FORMATS = ["txt","xlsx","tsv","csv"]
        cleaned_data = super(QueryMrna, self).clean()


        if not cleaned_data.get('matrix_mrna') and not cleaned_data.get('url'):
             self.add_error('matrix_mrna', 'One of these two fields is required (file or url)')
             self.add_error('url_mrna', 'One of these two fields is required (file or url)')
        if cleaned_data.get('matrix_mrna') and cleaned_data.get('url'):
            self.add_error('matrix_mrna', 'Choose either file or URL')
            self.add_error('url_mrna', 'Choose either file or URL')

        if cleaned_data.get('matrix_mrna'):
            fileType = cleaned_data.get('matrix_mrna').name.split(".")[1]
            if not fileType in ALLOWED_FORMATS:
                self.add_error('matrix_mrna', fileType+' is not an allowed matrix format. Please choose another file.')
        

        return cleaned_data


class QueryTrna(forms.Form):
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


    matrix_trna = forms.FileField(label='Upload input matrix', required=False)
    url_trna = forms.URLField(label='URL/link', required=False)
    jobID = forms.CharField(widget = forms.HiddenInput())
    typeJob = forms.CharField(widget = forms.HiddenInput(), initial = "miRNA")
    methods_trna = forms.MultipleChoiceField(choices=METHODS,label="Choose one or more normalization methods:",required=True)
    mirgenedb_trna = forms.ChoiceField(choices=METHODS,label="Choose short name from MiRGeneDB:",required=False)
    annotation_trna = forms.FileField(label="Annotation:", required=False)
    annotationURL_trna = forms.URLField(label='URL/link to annotation file:', required=False)
    minrc_trna = forms.IntegerField(label="Minimum RC",validators=[MinValueValidator(0)],initial=0,required=False)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.layout = Layout(
            create_card(
                TabHolder(
                        Tab(' Upload ',
                            HTML("""<br>"""),
                            Field('matrix_trna', css_class='form-control'),
                            css_id='upload_trna'
                        ),
                        Tab(' URL ',
                            HTML("""<br>"""),
                            Field('url_trna', css_class='form-control',placeholder="e.g.:https://url.com/matrix.txt"),
                            css_class="link",
                            css_id='url_trna'),
                        css_id='tab_upload_trna'
                ),
                title='Input Matrix',
                id="input_trna",
                show=True,
                modal="Choose_Input"
            ),
            HTML("""<br>"""),
            create_card(
                TabHolder(
                    Tab('Upload Annotation File',
                    HTML("""<br>"""),
                    Field('annotation_trna',css_class='form-control'),
                    css_id='annotation_upload_trna'
                    ),
                    Tab('URL of Annotation File',
                    HTML("""<br>"""),
                    Field('annotationURL_trna',css_class='form-control'),
                    css_id = 'annotation_url_trna'
                    ),
                    css_id="tab_annotation_trna"
                ),
                HTML('<a href="'+STATIC_URL+'testFiles/template.txt" download="template.txt"><button type="button" class="btn btn-info btn-sm float-right">Download sample template</button></a>'),
                title="Annotation",
                id="annotation_trna",
                show=True,
                modal="Choose_Annotation"
            ),
            HTML("""<br>"""),
            create_card(
                Field('minrc_trna',css_class='form-control'),
                title="Parameters",
                id="parameters_trna",
                show=False,
                modal="Parameters"
            ),
            HTML("""<br>"""),
            create_card(
                Field('methods_trna'),
                title= "Normalization Methods",
                id="normalization_trna",
                show=True,
                modal="Choose_Method"
            ),
            HTML("""<br>"""),
            'jobID',
            'typeJob',
            ButtonHolder(
                Submit('submit_trna', 'SUBMIT', css_class='btn btn-success btn-xl', onclick='check_form();')
            )
        )
    def clean(self):
        ALLOWED_FORMATS = ["txt","xlsx","tsv","csv"]
        cleaned_data = super(QueryTrna, self).clean()


        if not cleaned_data.get('matrix_trna') and not cleaned_data.get('url'):
             self.add_error('matrix_trna', 'One of these two fields is required (file or url)')
             self.add_error('url_trna', 'One of these two fields is required (file or url)')
        if cleaned_data.get('matrix_trna') and cleaned_data.get('url'):
            self.add_error('matrix_trna', 'Choose either file or URL')
            self.add_error('url_trna', 'Choose either file or URL')

        if cleaned_data.get('matrix_trna'):
            fileType = cleaned_data.get('matrix_trna').name.split(".")[1]
            if not fileType in ALLOWED_FORMATS:
                self.add_error('matrix_trna', fileType+' is not an allowed matrix format. Please choose another file.')
        

        return cleaned_data


class QueryOther(forms.Form):
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


    matrix_other = forms.FileField(label='Upload input matrix', required=False)
    url_other = forms.URLField(label='URL/link', required=False)
    jobID = forms.CharField(widget = forms.HiddenInput())
    typeJob = forms.CharField(widget = forms.HiddenInput(), initial = "other")
    methods_other = forms.MultipleChoiceField(choices=METHODS,label="Choose one or more normalization methods:",required=True)
    mirgenedb_other = forms.ChoiceField(choices=METHODS,label="Choose short name from MiRGeneDB:",required=False)
    annotation_other = forms.FileField(label="Annotation:", required=False)
    annotationURL_other = forms.URLField(label='URL/link to annotation file:', required=False)
    minrc_other = forms.IntegerField(label="Minimum RC",validators=[MinValueValidator(0)],initial=0,required=False)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.layout = Layout(
            create_card(
                TabHolder(
                        Tab(' Upload ',
                            HTML("""<br>"""),
                            Field('matrix_other', css_class='form-control'),
                            css_id='upload_other'
                        ),
                        Tab(' URL ',
                            HTML("""<br>"""),
                            Field('url_other', css_class='form-control',placeholder="e.g.:https://url.com/matrix.txt"),
                            css_class="link",
                            css_id='url_other'),
                        css_id='tab_upload_other'
                ),
                title='Input Matrix',
                id="input_other",
                show=True,
                modal="Choose_Input"
            ),
            HTML("""<br>"""),
            create_card(
                TabHolder(
                    Tab('Upload Annotation File',
                    HTML("""<br>"""),
                    Field('annotation_other',css_class='form-control'),
                    css_id='annotation_upload_other'
                    ),
                    Tab('URL of Annotation File',
                    HTML("""<br>"""),
                    Field('annotationURL_other',css_class='form-control'),
                    css_id = 'annotation_url_other'
                    ),
                    css_id="tab_annotation_other"
                ),
                HTML('<a href="'+STATIC_URL+'testFiles/template.txt" download="template.txt"><button type="button" class="btn btn-info btn-sm float-right">Download sample template</button></a>'),
                title="Annotation",
                id="annotation_other",
                show=True,
                modal="Choose_Annotation"
            ),
            HTML("""<br>"""),
            create_card(
                Field('minrc_other',css_class='form-control'),
                title="Parameters",
                id="parameters_other",
                show=False,
                modal="Parameters"
            ),
            HTML("""<br>"""),
            create_card(
                Field('methods_other'),
                title= "Normalization Methods",
                id="normalization_other",
                show=True,
                modal="Choose_Method"
            ),
            HTML("""<br>"""),
            'jobID',
            'typeJob',
            ButtonHolder(
                Submit('submit_other', 'SUBMIT', css_class='btn btn-success btn-xl', onclick='check_form();')
            )
        )
    def clean(self):
        ALLOWED_FORMATS = ["txt","xlsx","tsv","csv"]
        cleaned_data = super(QueryOther, self).clean()


        if not cleaned_data.get('matrix_other') and not cleaned_data.get('url'):
             self.add_error('matrix_other', 'One of these two fields is required (file or url)')
             self.add_error('url_other', 'One of these two fields is required (file or url)')
        if cleaned_data.get('matrix_other') and cleaned_data.get('url'):
            self.add_error('matrix_other', 'Choose either file or URL')
            self.add_error('url_other', 'Choose either file or URL')

        if cleaned_data.get('matrix_other'):
            fileType = cleaned_data.get('matrix_other').name.split(".")[1]
            if not fileType in ALLOWED_FORMATS:
                self.add_error('matrix_other', fileType+' is not an allowed matrix format. Please choose another file.')
        

        return cleaned_data

