#!/usr/bin/python2
  
## /home/jhy/spatial_transcriptomics/220905/code/step2_Seurat_preprocessing.py /home/jhy/spatial_transcriptomics/220905/result/spaceranger/CPS_OV10RTOV4-1/CPS_OV10RTOV4-1/outs/ /home/jhy/spatial_transcriptomics/220905/result/seurat/ OV_Spatial

## /home/jhy/spatial_transcriptomics/220905/code/step2_Seurat_preprocessing.py /home/jhy/spatial_transcriptomics/220905/result/spaceranger/ /home/jhy/spatial_transcriptomics/220905/result/seurat/

import sys, os
from subprocess import Popen, PIPE, STDOUT
import time, commands
import parameter as pr

inputdir = sys.argv[1]
outputdir = sys.argv[2]

samples = os.listdir(inputdir)
samples = filter(lambda x: not x.startswith('slurm'), samples)
#print samples

for sample in samples:
    indir = inputdir+sample+'/'+sample+'/outs/'
    otdir = outputdir+sample+'/'
    #pr.make_dir(otdir)
    project = sample
    jobname = 'ST'
    commandline = 'sbatch -J ' + jobname + ' -D ' + otdir + '/ -p gar --wrap="Rscript /home/jhy/spatial_transcriptomics/220905/code/Seurat_preprocessing.R ' + indir + ' ' + otdir + ' ' + project + '"'
    pr.inShell(commandline)


