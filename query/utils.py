from crispy_forms.layout import HTML, Div


def create_card(*fields, **kwargs):
    title = kwargs['title']
    id_card = kwargs['id']
    show = kwargs['show']
    modal = kwargs['modal']


    card =  HTML("<div class='card'>")
    panel_title_1 = HTML("<div class='card-header card-query'>")
    panel_title_2 = HTML("<div class='row'><div class='col-9'><button class='btn btn-link btn-block text-left' type='button' data-toggle='collapse' data-target='#"+id_card+"' aria-expanded='true' aria-controls='"+id_card+"'><h5>"+title+"</h5></button></div>")
    if modal!=None:
        panel_title_3 = HTML("<div class='col-3'><div class='float-right'>"+render_modal(modal)+"</div></div></div></div>")
    else:
        panel_title_3 = HTML("<div class='col-3'><div class='float-right'></div></div></div></div>")
      
    panel_title = Div(panel_title_1,panel_title_2,panel_title_3)
    if show:
        card_body=HTML("<div id='"+id_card+"' class='collapse show'><div class='card-body'>")
    else:
        card_body=HTML("<div id='"+id_card+"' class='collapse'><div class='card-body'>")
    form_group = Div(*fields, css_class='form-group')
    body = Div(form_group)
    end = HTML("</div></div></div>")
    panel = Div(card, panel_title,card_body,body,end)
    return panel

def render_modal(modal):
    return '<a  data-toggle="modal" data-target="#'+modal+'" class="btn btn-link"><i class="fa fa-question-circle"> </i></a>'

def render_modal_annotation(modal):
    return '<a  data-toggle="modal" data-target="#'+modal+'" class="btn btn-link"><button type="button" class="btn btn-secondary btn-sm">Upload Annotation</button></a>'

def create_collapsable_div(*fields, **kwargs):
    title = kwargs['title']
    c_id = kwargs['c_id']
    extra = ''
    if 'extra_title' in kwargs:
        extra = kwargs['extra_title']
    panel_title = HTML('<h4 class="panel-title"><a data-toggle="collapse" href="#'+c_id+'">'+title+'</a>'+extra+'</h4>')
    panel_heading = Div(panel_title, css_class='panel-heading')
    form_group = Div(*fields, css_class='form-group')
    body = Div(form_group, css_class='panel-body')
    if 'open' in kwargs and kwargs['open']:
        collapse_in = "panel-collapse collapse in"
    else:
        collapse_in = "panel-collapse collapse"
    panel_collapse = Div(body, css_class=collapse_in, css_id=c_id)
    panel = Div(panel_heading, panel_collapse, css_class="panel panel-default")
    return panel