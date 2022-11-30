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

# Make normalization and plots

methods = config['methods']

for method in methods:
    norm(df,method,jobDir)
    createplots(df,method,jobDir)
    

#Check and create results.txt

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)

sys.exit(0)