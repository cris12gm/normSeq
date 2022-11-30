from crispy_forms.layout import HTML, Div


def create_card(*fields, **kwargs):
    title = kwargs['title']

    card =  HTML("<div class='card'>")
    panel_title = HTML("<div class='card-header'><h5>"+title+"</h5></div>")
    card_body=HTML("<div class='card-body'>")
    form_group = Div(*fields, css_class='form-group')
    body = Div(form_group)
    end = HTML("</div></div>")
    panel = Div(card, panel_title,card_body,body,end)
    return panel

def render_modal(modal):
    return '<a  data-toggle="modal" data-target="#'+modal+'" class="btn btn-link"><i class="fa fa-question-circle"> </i></a>'

def create_tabs(*fields,**kwargs):
    tabs = kwargs['tabs']

    form_group = Div(*fields, css_class='form-group')
    header = HTML('<ul class="nav nav-tabs" id="myTab" role="tablist">')
    endHeader = HTML('</ul>')
    headerTabContent = HTML('<div class="tab-content" id="myTabContent">')
    endHeaderTabContent = HTML('</div>')
    nTabs = 0
    for tab in tabs:
        if nTabs==0:
            contentHeader = HTML('<li class="nav-item"><a class="nav-link active" id="'+tab+'-tab" data-toggle="tab" href="#'+tab+'" role="tab" aria-controls="'+tab+'" aria-selected="true">'+tab+'</a></li>')
            contentTab = HTML('<div class="tab-pane fade show active" id="'+tab+'" role="tabpanel" aria-labelledby="'+tab+'-tab">...</div>')
        else:
            contentHeader = contentHeader + HTML('<li class="nav-item"><a class="nav-link active" id="'+tab+'-tab" data-toggle="tab" href="#'+tab+'" role="tab" aria-controls="'+tab+'" aria-selected="false">'+tab+'</a></li>')
            contentTab = contentTab + HTML('<div class="tab-pane fade show active" id="'+tab+'" role="tabpanel" aria-labelledby="'+tab+'-tab">...</div>')

    panel = Div(header,contentHeader,endHeader,headerTabContent,contentTab,endHeaderTabContent)
    return panel
