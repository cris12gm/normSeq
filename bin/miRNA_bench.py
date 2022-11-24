import os,sys
from plots import heatmap
from miRNAnorm import readTxt

jobDir = sys.argv[1]
# Make normalization

infile = os.path.join(jobDir,"matrix.txt")
outfile = os.path.join(jobDir,"matrix_RPM.txt")

#print(infile)
readTxt(infile,outfile)

# Create heatmap

infile = os.path.join(jobDir,"matrix.txt")
outfile = os.path.join(jobDir,"heatmap.html")

heatmap(infile,outfile)


#Check and create results.txt

resultsFile = os.path.join(jobDir,"results.txt")
os.system("touch "+resultsFile)