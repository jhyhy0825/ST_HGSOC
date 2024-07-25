#!/usr/bin/python2
  
## /home/jhy/spatial_transcriptomics/220905/code/Giotto/Giotto4.py /home/jhy/spatial_transcriptomics/220905/result/seurat/

import sys, os
from subprocess import Popen, PIPE, STDOUT
import time, commands
import parameter as pr

inputdir = sys.argv[1]

samples = os.listdir(inputdir)
samples = filter(lambda x: not x.startswith('nothing'), samples)
samples = filter(lambda x: not x.startswith('gene_list'), samples)

for sample in samples:
    job = 'ST'
    sampledir = inputdir+'/'+sample+'/'
    outputdir = inputdir+'/'+sample+'/Giotto_GSE180661_matrix/'
    pr.make_dir(outputdir)
    commandline = 'sbatch -J ' + job + ' -D ' + outputdir + '/ -p gar --wrap="Rscript /home/jhy/spatial_transcriptomics/220905/code/Giotto/Giotto4.R /home/jhy/single_cell_lung_mouse/HumanValidation/rawdata/GSE180661/Adnexa_count_matrix.txt.gz /home/jhy/single_cell_lung_mouse/HumanValidation/rawdata/GSE180661/Adnexa_cell_type.txt ' + outputdir + ' ' +  sampledir + '"'
    pr.inShell(commandline)
    while (commands.getoutput('squeue -n' + job).count(job)>2): time.sleep(10)
    time.sleep(1)
