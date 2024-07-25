#!/usr/bin/python2
  
## /home/jhy/spatial_transcriptomics/220905/code/step3_Seurat_dimensional_reduction_combine.py /home/jhy/spatial_transcriptomics/220905/result/seurat/

import sys, os
from subprocess import Popen, PIPE, STDOUT
import time, commands
import parameter as pr

inputdir = sys.argv[1]

samples = os.listdir(inputdir)
samples = filter(lambda x: not x.startswith('nothing'), samples)
samples = ['combine']

for sample in samples:
    otdir = inputdir+sample+'/'
    jobname = 'ST'
    commandline = 'sbatch -J ' + jobname + ' -D ' + otdir + '/ -p gar --wrap="Rscript /home/jhy/spatial_transcriptomics/220905/code/Seurat_dimensional_reduction.R ' + otdir + '"'
    pr.inShell(commandline)


