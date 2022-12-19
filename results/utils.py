import math

def creategrid(content):
    init = '<div class="grid">'
    end = '</div>'

    divobject = init+content+end
    return divobject

def createrow(content):
    init = '<div class="row">'
    end = '</div>'

    divobject = init+content+end

    return divobject

def createcol(content,typeCol):
    init = '<div class="'+typeCol+'">'
    end = '</div>'

    divobject = init+content+end

    return divobject

def createCardImage(content):
    init = '<div class="card" style="width: 100%;"><div class="row"><div class="col-10 my-auto mx-auto">'
    imageCode = '<img class="img-fluid" src="'+content+'">'
    end = '</div></div></div>'

    divobject = init+imageCode+end
    return divobject

def createGridPlot(dictData):
    numRows = math.ceil(len(dictData) / 3)
    for i in range(0,numRows,1):
        thisRow = ''
        for i in range(0,3,1):
            thisRow = createcol(createCardImage(dictData[0]))