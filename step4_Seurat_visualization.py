#!/usr/bin/python2
  
## /home/jhy/spatial_transcriptomics/220905/code/step4_Seurat_visualization.py /home/jhy/spatial_transcriptomics/220905/result/seurat/

import sys, os
from subprocess import Popen, PIPE, STDOUT
import time, commands
import parameter as pr

inputdir = sys.argv[1]

samples = os.listdir(inputdir)
samples = filter(lambda x: not x.startswith('nothing'), samples)
#samples = ['CPS_OV34RtOV1']

for sample in samples:
    otdir = inputdir+sample+'/'
    jobname = 'ST'
    commandline = 'sbatch -J ' + jobname + ' -D ' + otdir + '/ -p gar --wrap="Rscript /home/jhy/spatial_transcriptomics/220905/code/Seurat_visualization.R ' + otdir + '"'
    pr.inShell(commandline)


