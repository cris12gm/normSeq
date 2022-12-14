#C:\Users\Cris\Dropbox\TRABAJO\NL\normSeq\Website\myvenv\ python

import os,sys
import json
from plots import createplots
from miRNAnorm import processInput,norm


#Error = False

configFile = open(sys.argv[1])
config = json.load(configFile)

jobDir = config['jobDir']

#Get df
df = processInput(os.path.join(jobDir,"matrix.txt"))
infile = os.path.join(jobDir,"matrix.txt")
annotation = os.path.join(jobDir,"annotation.txt")

# Make normalization and plots

methods = config['methods']

for method in methods:

    outdf,normfile = norm(infile,df,method,jobDir)
    createplots(normfile,outdf,method,jobDir,annotation)
    

#Check and create results.txt

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)

sys.exit(0)