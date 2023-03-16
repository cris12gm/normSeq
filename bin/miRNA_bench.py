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
import itertools

from subprocess import Popen
from multiprocessing import Pool

configFile = open(sys.argv[1])
config = json.load(configFile)

jobDir = config['jobDir']

#Error, status and log Files
logFile = os.path.join(jobDir,"Log.txt")
errorFile = os.path.join(jobDir,"Error.txt")
statusFile = os.path.join(jobDir,"status.txt")

log = open(logFile, 'a')
error = open(errorFile, 'a')
status = open(statusFile, 'a')

#Get df
df = processInput(os.path.join(jobDir,"matrix.txt"))

#Create directories

##Normalized
if not os.path.exists(os.path.join(jobDir,"normalized")):
    os.mkdir(os.path.join(jobDir,"normalized"))

#Save df in normalized
outfile_NN = os.path.join(jobDir,"normalized","matrix_NN.txt") 
df.to_csv(outfile_NN,sep="\t")

if os.path.isfile(outfile_NN):
    log.write("0. No normalized file saved\n")
else:
    error.write("Normalized file couldn't be safe\n")

#Define files for the following steps
infile = os.path.join(jobDir,"normalized","matrix_NN.txt")

annotation = os.path.join(jobDir,"annotation.txt")
try:
    annotation_df=processAnnotation(annotation)
    os.system("mv "+annotation+" "+jobDir+"/annotation_old.txt")
    annotation_df.to_csv(annotation,sep="\t")
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

##################################################################
################## DIFFERENTIAL EXPRESSION #######################
##################################################################

if config['diffExpr']=="True":
    FDR = config['pval']
    min_t = str(0)
    methodDE = "TMM" #Change in input
    if not os.path.exists(os.path.join(jobDir,"DE")):
        os.mkdir(os.path.join(jobDir,"DE"))
    combinations = createGroupFile(annotation_df,jobDir)

    output_de = de_R(infile,annotation,combinations,methodDE,FDR,min_t,jobDir,error,log,status)
    output_de = ttest(df,combinations,annotation_df,FDR,output_de,jobDir,error,log,status)
    consensus(output_de,jobDir)

##################################################################
####################### NORMALIZATION ############################
##################################################################

#Methods from R code -> Parallelize with os

methods_r = ["UQ","TMM","RLE","DESEQ","QN","RUV"]
cmds_r = []
r_files = []

#Dict with all the normalized datasets
normalized = {}

log.write("3. Normalization\n")
status.write("<p>3. Normalization</p>")
status.flush()

for method in methods:
    #Normalization of non R + save R norms in cmds_r
    if method in methods_r:
        cmds_r,outfile = norm_r(infile,method,jobDir,cmds_r)
        r_files.append([outfile,method])
    else:
        outdf,normfile = norm(infile,df,method,jobDir)
        normalized[method] = [outdf,normfile]

#Launch the normalization of Rs
procs = [ Popen(i,shell=True) for i in cmds_r ]
for p in procs:
    p.wait()

#Create the matrix from the Rs
for file,method in r_files:
    normfile = file
    outdf = processInput(normfile)
    normalized[method] = [outdf,normfile]

##################################################################
##################### GRAPHS DIR CREATION ########################
##################################################################

graphsDir = os.path.join(jobDir,"graphs")
summaryDir = os.path.join(jobDir,"graphs","summary")
if not os.path.exists(graphsDir):
    os.mkdir(graphsDir)
if not os.path.exists(summaryDir):
    os.mkdir(summaryDir)

##################################################################
####################### GROUPS CREATION ##########################
##################################################################

groups = annotation_df["group"].values.tolist()
diffGroups = list(set(groups))
combinations = []
output = open(os.path.join(summaryDir,"groups.txt"),'w')
for subset in itertools.combinations(diffGroups, 2):
    combinations.append([subset[0],subset[1]])
    element = subset[0]+"-"+subset[1]
    output.write(element+"\n")
output.close()

##################################################################
###################### #INFORMATION GAIN #########################
##################################################################
log.write("4. Information Gain analysis\n")
status.write("<p>4. Information Gain analysis</p>")
status.flush()

for combination in combinations:
    comb = combination[0]+"-"+combination[1]
    information_gain = {}
    for method in methods:
        normdf = normalized[method][0]
        info = calculate_infoGain(normdf,annotation_df,combination)
        information_gain[method] = info
    outfileImage = os.path.join(jobDir,"graphs","summary","infoGain_"+combination[0]+"-"+combination[1]+".png")
    outfile = os.path.join(jobDir,"graphs","summary","infoGain_"+combination[0]+"-"+combination[1]+".html")
    title = combination[0]+"-"+combination[1]
    plotInfo(information_gain,outfileImage,outfile,title)

##################################################################
####################### SUMMARY + PLOTS ##########################
##################################################################

cmd_plots = []
status.write("<p>5. Visualizations</p>")
status.flush()

for method in normalized:
    normdf = normalized[method][0]
    normfile = normalized[method][1]

#Summary
    createsummary(normfile,normdf,method,jobDir,annotation_df,combinations)
    
    log.write("5. Visualization\n")

#Plots
    cmd_plots = createplots(normfile,normdf,method,jobDir,annotation,annotation_df,cmd_plots)

#Launch the Rs from plots

procs = [ Popen(i,shell=True) for i in cmd_plots ]
for p in procs:
    p.wait()

##################################################################
####################### FINAL CHECK ##############################
##################################################################

#Check and create results.txt
log.close()
error.close()
status.close()

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)

sys.exit(0)