import os,sys
import pandas as pd
import numpy as np
import subprocess

from config import R_PATH
from config import R_SCRIPTS_PATH

def combat(infile,annotation):
    subprocess.call (R_PATH+" --vanilla "+R_SCRIPTS_PATH+"combat.R "+infile+" "+annotation,shell=True)
    sys.exit(1)