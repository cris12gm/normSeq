#C:\Users\Cris\Dropbox\TRABAJO\NL\normSeq\Website\myvenv\ python

import os,sys
import json
from plots import createplots
from miRNAnorm import processInput,norm
from summary import createsummary
from correctBatch import combat,plotsBatch
from de import edgeR,processAnnotation
from infoGain import calculate_infoGain,plotInfo


#Error = False

configFile = open(sys.argv[1])
config = json.load(configFile)

jobDir = config['jobDir']

#Get df
df = processInput(os.path.join(jobDir,"matrix.txt"))
infile = os.path.join(jobDir,"matrix.txt")
annotation = os.path.join(jobDir,"annotation.txt")

#Correct batch effect
if config['batchEffect'] == 'True':
    batchAnnotation = os.path.join(jobDir,"batchFile.txt")
    outfile = os.path.join(jobDir,"matrix_corrected.txt")
    combat(infile,batchAnnotation,outfile)
    oldMatrix = os.path.join(jobDir,"matrix_old.txt")
    os.rename(infile, oldMatrix)
    os.rename(outfile, infile)

    dfCorrected = processInput(os.path.join(jobDir,"matrix.txt"))
    dfOld = processInput(os.path.join(jobDir,"matrix_old.txt"))

    plotsBatch(dfCorrected,dfOld,annotation,jobDir)


# Make normalization and plots
methods = config['methods']

infoGain = {}
for method in methods:

    #if method=='UQ' or method=='TMM' or method=='RLE':
    #    FDR = config['pval']
    #    edgeR(infile,method,annotation,FDR,jobDir)
    outdf,normfile = norm(infile,df,method,jobDir)

    criterion='entropy'
    info_method = calculate_infoGain(normfile,annotation,criterion)
    infoGain[method] = info_method
    createsummary(normfile,outdf,method,jobDir,annotation)
    createplots(normfile,outdf,method,jobDir,annotation)

outfileImage = os.path.join(jobDir,"graphs","summary","infoGain.png")
outfile = os.path.join(jobDir,"graphs","summary","infoGain.html")
plotInfo(infoGain,outfileImage,outfile)


#Check and create results.txt

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)

sys.exit(0)