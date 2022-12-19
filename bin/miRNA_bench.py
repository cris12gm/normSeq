#C:\Users\Cris\Dropbox\TRABAJO\NL\normSeq\Website\myvenv\ python

import os,sys
import json
from plots import createplots
from miRNAnorm import processInput,norm
from summary import createsummary
from correctBatch import combat,plotsBatch


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

for method in methods:

    outdf,normfile = norm(infile,df,method,jobDir)
    createsummary(normfile,outdf,method,jobDir,annotation)
    createplots(normfile,outdf,method,jobDir,annotation)

    

#Check and create results.txt

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)

sys.exit(0)