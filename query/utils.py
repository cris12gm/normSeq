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
