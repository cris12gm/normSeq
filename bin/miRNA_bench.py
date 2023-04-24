#C:\Users\Cris\Dropbox\TRABAJO\NL\normSeq\Website\myvenv\ python

import os,sys
import json
from plots import createplots
from miRNAnorm import processInput,norm,norm_r,processAnnotation,processInputInit,rleplot
from summary import createsummary
from correctBatch import combat,plotsBatch
from de import de_R,createGroupFile,ttest,consensus,plotDE
from infoGain import calculate_infoGain_group,calculate_infoGain_pairwise,plotInfo
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
status = open(statusFile, 'a')

##################################################################
######################### ANNOTATION #############################
##################################################################

annotation = os.path.join(jobDir,"annotation.txt")
try:
    annotation_df,samples,replicates=processAnnotation(annotation)
    os.system("mv "+annotation+" "+jobDir+"/annotation_old.txt")
    annotation_df.to_csv(annotation,sep="\t")
    log.write("1. Annotation processed\n")
    status.write("<p>1. Input matrix and annotation has been processed</p>")
    status.flush()
except:
    error = open(errorFile, 'a')
    error.write("<p>There is an error in the annotation file format. Please, check the format guidelines <a href='https://arn.ugr.es/normseq_doc/annotation/'>here</a></p>")
    error.close()
    sys.exit(0)

##################################################################
######################### INPUT DATA #############################
##################################################################

#Get df
try:
    df,dfOriginal = processInputInit(os.path.join(jobDir,"matrix.txt"),samples,config['minRC'],annotation_df)
except:
    error = open(errorFile, 'a')
    error.write("<p>Normalized file couldn't be safe. Please, check the format guidelines <a href='https://arn.ugr.es/normseq_doc/annotation/'>here</a></p>")
    error.close()
    sys.exit(0)

#Create directories
##Normalized
if not os.path.exists(os.path.join(jobDir,"normalized")):
    os.mkdir(os.path.join(jobDir,"normalized"))

#Save df in normalized
outfile_NN = os.path.join(jobDir,"normalized","matrix_NN.txt") 
df.to_csv(outfile_NN,sep="\t")
outfile_NN_Original = os.path.join(jobDir,"normalized","matrix_NN_toDE.txt") 
dfOriginal.to_csv(outfile_NN_Original,sep="\t")

if os.path.isfile(outfile_NN):
    log.write("0. No normalized file saved\n")
else:
    error = open(errorFile, 'a')
    error.write("<p>Normalized file couldn't be safe. Please, check the format guidelines <a href='https://arn.ugr.es/normseq_doc/annotation/'>https://arn.ugr.es/normseq_doc/annotation/>here</a></p>")
    error.close()
    sys.exit(0)

#Define files for the following steps
infile = os.path.join(jobDir,"normalized","matrix_NN.txt")

##################################################################
######################## BATCH EFFECT ############################
##################################################################

if config['batchEffect'] == 'True':
    batchAnnotation = os.path.join(jobDir,"batchFile.txt")
    outfile = os.path.join(jobDir,"matrix_corrected.txt")
    combat(infile,batchAnnotation,annotation,outfile)
    matrixFile = os.path.join(jobDir,"matrix.txt")
    oldMatrix = os.path.join(jobDir,"matrix_old.txt")
    os.rename(matrixFile, oldMatrix)
    os.rename(outfile, infile)
    os.system("cp "+infile+" "+matrixFile)
    os.system("cp "+matrixFile+" "+infile)
    try:
        dfCorrected = processInput(os.path.join(jobDir,"matrix.txt"),annotation_df)
        outfile_NN_Original = os.path.join(jobDir,"normalized","matrix_NN_toDE.txt") 
        dfCorrected.to_csv(outfile_NN_Original,sep="\t")
        dfOld = processInput(os.path.join(jobDir,"matrix_old.txt"),annotation_df)
    except:
        error = open(errorFile, 'a')
        error.write("<p>Batch effect correction was not possible, please check the input files. Please, check the format guidelines <a href='https://arn.ugr.es/normseq_doc/annotation/'>https://arn.ugr.es/normseq_doc/annotation/>here</a></p>")
        error.close()
        sys.exit(0)
    
    batch_df,samples,discard = processAnnotation(batchAnnotation)
    plotsBatch(dfCorrected,dfOld,batch_df,jobDir)
    df = dfCorrected

    os.system("rm "+outfile_NN)
    df.to_csv(outfile_NN,sep="\t")

##################################################################
################## DIFFERENTIAL EXPRESSION #######################
##################################################################

if config['diffExpr']=="True":
    FDR = config['pval']
    min_t = str(0)
    methodDE = "TMM" #Change in input
    methodologyEdgeR = config['edgeR']
    if not os.path.exists(os.path.join(jobDir,"DE")):
        os.mkdir(os.path.join(jobDir,"DE"))
    combinations = createGroupFile(annotation_df,jobDir)

    output_de = de_R(outfile_NN_Original,annotation,combinations,methodDE,FDR,min_t,methodologyEdgeR,jobDir,log,status)
    output_de = ttest(df,combinations,annotation_df,FDR,output_de,jobDir,log,status)
    consensus(output_de,df,annotation_df,jobDir)
    plotDE(df,output_de,annotation_df,jobDir)

##################################################################
####################### NORMALIZATION ############################
##################################################################
# Make normalization and plots
methods = config['methods']
infoGain = {}

#Methods from R code -> Parallelize with os

methods_r = ["TMM","RLE","DESEQ","QN","RUV"]
cmds_r = []
r_files = []

#Dict with all the normalized datasets
normalized = {}

log.write("3. Normalization\n")
status.write("<p>3. Normalization</p>")
status.flush()

for method in methods:
    if method=="RUV" and replicates==False:
        pass
    #Normalization of non R + save R norms in cmds_r
    if method in methods_r:
        cmds_r,outfile = norm_r(infile,method,jobDir,annotation,cmds_r)
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
    if method=='RUV' and replicates==False:
        pass
    outdf = processInput(normfile,annotation_df)
    normalized[method] = [outdf,normfile]


#RLE plots

cmds_rle = []
for method in methods:
    cmd_method = rleplot(normalized[method][1],annotation,jobDir,method)
    cmds_rle.append(cmd_method)

procs = [ Popen(i,shell=True) for i in cmds_rle ]
for p in procs:
    p.wait()

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
output = open(os.path.join(summaryDir,"groups.txt"),'w')
for group in diffGroups:
    output.write(group+"\n")
output.close()


output = open(os.path.join(summaryDir,"combinations.txt"),'w')
combinations = []
for subset in itertools.combinations(diffGroups, 2):
    combinations.append([subset[0],subset[1]])
    output.write(subset[0]+"-"+subset[1]+"\n")

output.close()


##################################################################
####################### INFORMATION GAIN #########################
##################################################################
log.write("4. Information Gain analysis\n")
status.write("<p>4. Information Gain analysis</p>")
status.flush()


#Per Group
for group in diffGroups:
    information_gain = {}
    for method in methods:
        normdf = normalized[method][0]
        info = calculate_infoGain_group(normdf,annotation_df,group)
        information_gain[method] = info
    # infoDf = pd.DataFrame(information_gain)
    infoDf = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in information_gain.items() ]))
    infoDf['name'] = list(normdf.index)
    infoDf = infoDf.set_index('name')

    # # #Get Top 10
    # # outfileImage = os.path.join(jobDir,"graphs","summary","infoGainTop10_"+group+".txt")
    # # outfile = os.path.join(jobDir,"graphs","summary","infoGainTop10_"+group+".html")
    # # top10Info(infoDf,normalized,methods,annotation_df,outfileImage,outfile)
    # # sys.exit(1)
    outfileInfo = os.path.join(jobDir,"graphs","summary","infoGain_"+group+".txt")
    infoDf.to_csv(outfileInfo,sep="\t")
    outfileImage = os.path.join(jobDir,"graphs","summary","infoGain_"+group+".png")
    outfile = os.path.join(jobDir,"graphs","summary","infoGain_"+group+".html")
    title = group
    titleaxis = "Information Gain per "+config["typeJob"]
    plotInfo(information_gain,outfileImage,outfile,title,titleaxis)

#Pairwise

for combination in combinations:
    information_gain_groups = {}
    for method in methods:
        normdf = normalized[method][0]
        info = calculate_infoGain_pairwise(normdf,annotation_df,combination)
        information_gain_groups[method] = info
    outfileImage = os.path.join(jobDir,"graphs","summary","infoGain_"+combination[0]+"-"+combination[1]+".png")
    outfile = os.path.join(jobDir,"graphs","summary","infoGain_"+combination[0]+"-"+combination[1]+".html")
    title = combination[0]+"-"+combination[1]
    titleaxis = "Information Gain per "+config["typeJob"]
    plotInfo(information_gain_groups,outfileImage,outfile,title,titleaxis)

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
status.close()

os.system("rm "+jobDir+"/annotation_old.txt")

zipFile = os.path.join(jobDir,"results_"+config["jobID"]+".zip")
zipCmd = "cd "+jobDir+" && zip -r "+zipFile+" ./*"

os.system(zipCmd)

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)


sys.exit(0)