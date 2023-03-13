#C:\Users\Cris\Dropbox\TRABAJO\NL\normSeq\Website\myvenv\ python

import os,sys
import json
from plots import createplots
from miRNAnorm import processInput,norm,norm_r,processAnnotation
from summary import createsummary
from correctBatch import combat,plotsBatch
from de import de_R,createGroupFile,ttest,consensus
from infoGain import calculate_infoGain,plotInfo
import pandas as pd

from subprocess import Popen


#Error = False

configFile = open(sys.argv[1])
config = json.load(configFile)

jobDir = config['jobDir']

#Error, status and log Files
logFile = os.path.join(jobDir,"Log.txt")
errorFile = os.path.join(jobDir,"Error.txt")
statusFile = os.path.join(jobDir,"status.txt")

log = open(logFile, 'w')
error = open(errorFile, 'w')
status = open(statusFile, 'w')

#Get df
df = processInput(os.path.join(jobDir,"matrix.txt"))

#Create directories

##Normalized
if not os.path.exists(os.path.join(jobDir,"normalized")):
    os.mkdir(os.path.join(jobDir,"normalized"))

#Save df in normalized
outfile_NN = os.path.join(jobDir,"normalized","matrix_NN.txt") 

if os.path.isfile(outfile_NN):
    df.to_csv(outfile_NN,sep="\t")
    log.write("0. No normalized file saved\n")
else:
    error.write("Normalized file couldn't be safe\n")

#Define files for the following steps
infile = os.path.join(jobDir,"normalized","matrix_NN.txt")

annotation = os.path.join(jobDir,"annotation.txt")
try:
    annotation_df=processAnnotation(annotation)
    log.write("1. Annotation processed\n")
    status.write("<p>1. Input matrix and annotation has been processed</p>")
    status.flush()
except:
    error.write("There is an error in the annotation file format. Please, check the format guidelines <a href='https://arn.ugr.es/normseq_doc/annotation/'>https://arn.ugr.es/normseq_doc/annotation/>here</a>")
    status.write("<p>There is an error in the annotation file format. Please, check the format guidelines <a href='https://arn.ugr.es/normseq_doc/annotation/'>https://arn.ugr.es/normseq_doc/annotation/>here</a></p>")
    status.flush()
    sys.exit(0)

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
    df = dfCorrected


# Make normalization and plots
methods = config['methods']
infoGain = {}

#Differential Expression

if config['diffExpr']=="True":
    FDR = config['pval']
    min_t = str(0)
    method = "TMM" #Change in input
    if not os.path.exists(os.path.join(jobDir,"DE")):
        os.mkdir(os.path.join(jobDir,"DE"))
    combinations = createGroupFile(annotation_df,jobDir)

    output_de = de_R(infile,annotation,combinations,method,FDR,min_t,jobDir,error,log,status)
    output_de = ttest(df,combinations,annotation_df,FDR,output_de,jobDir,error,log,status)
    consensus(output_de,jobDir)

cmds_r = []

#Normalization
methods_r = ["UQ","TMM","RLE","DESEQ","QN","RUV"]

normalized = {}
r_files = []

log.write("4. Normalization\n")
status.write("<p>4. Normalization</p>")
status.flush()

for method in methods:
    #Normalization
    if method in methods_r:
        cmds_r,outfile = norm_r(infile,method,jobDir,cmds_r)
        r_files.append([outfile,method])
    else:
        outdf,normfile = norm(infile,df,method,jobDir)
        normalized[method] = [outdf,normfile]

#Launch the Rs
procs = [ Popen(i,shell=True) for i in cmds_r ]
for p in procs:
    p.wait()

#Create the matrix from the Rs
for file,method in r_files:
    normfile = file
    outdf = processInput(normfile)
    normalized[method] = [outdf,normfile]

cmd_plots = []
for method in normalized:
    outdf = normalized[method][0]
    normfile = normalized[method][1]
#Information Gain
    log.write("5. Information Gain analysis\n")
    status.write("<p>5. Information Gain analysis</p>")
    criterion=config["infoGain"]
    info_method = calculate_infoGain(normfile,annotation,criterion)
    infoGain[method] = info_method

#Summary
    createsummary(normfile,outdf,method,jobDir,annotation_df)

#Plots
    log.write("6. Visualization\n")
    status.write("<p>6. Visualization</p>")
    status.flush()
    cmd_plots = createplots(normfile,outdf,method,jobDir,annotation,annotation_df,cmd_plots)
#Launch the Rs from plots
procs = [ Popen(i,shell=True) for i in cmd_plots ]
for p in procs:
    p.wait()


#Plot Info Gain
outfileImage = os.path.join(jobDir,"graphs","summary","infoGain.png")
outfile = os.path.join(jobDir,"graphs","summary","infoGain.html")
plotInfo(infoGain,outfileImage,outfile)


#Check and create results.txt
log.close()
error.close()

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)

sys.exit(0)