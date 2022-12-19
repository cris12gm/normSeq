import os,sys
import pandas as pd
import numpy as np
import subprocess
from plots import pca
from config import R_PATH
from config import R_SCRIPTS_PATH

def combat(infile,annotation,outfile):
    subprocess.call (R_PATH+" --vanilla "+R_SCRIPTS_PATH+"combat.R "+infile+" "+annotation+" "+outfile,shell=True)

def plotsBatch(dfCorrected,dfOld,annotation,jobDir):
    
    outDir = os.path.join(jobDir,"graphs","batchEffect")
    if not os.path.exists(os.path.join(jobDir,"graphs")):
        os.mkdir(os.path.join(jobDir,"graphs"))
    
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    outfile = os.path.join(outDir,"pca_corrected.html")
    outfileImage = os.path.join(outDir,"pca_corrected.png")
    pca(dfCorrected,outfile,outfileImage)
    outfile = os.path.join(outDir,"pca_old.html")
    outfileImage = os.path.join(outDir,"pca_old.png")
    pca(dfOld,outfile,outfileImage)