#!/usr/bin/python2
  
## /home/jhy/spatial_transcriptomics/220905/code/step1_spaceranger_count.py /home/jhy/spatial_transcriptomics/220905/

import sys, os, time, commands
from subprocess import Popen, PIPE, STDOUT
import parameter as pr

inputdir = sys.argv[1]
rawdir = inputdir+'/HN00174122/HN00174122_10X_RawData_Outs/'
refdir = inputdir+'/spaceranger/ref/refdata-gex-GRCh38-2020-A/'
resultdir = inputdir+'/result/spaceranger/'
pr.make_dir(resultdir)
#samplename = 'HN00174122'
imagedir = inputdir+'/HN00174122/HN00174122_Image/NewName/'
names = imagedir+'samplename_and_filename.txt'
with open(names, 'r') as file:
    lines = file.readlines()
name = dict()
for line in lines:
    sp = line.split('\t')
    name[sp[0]] = sp[1]

sampledirs = os.listdir(rawdir)
sampledirs = filter(lambda x: not x.startswith('RawData'), sampledirs)
for sampledir in sampledirs:
    fastq_path = rawdir+sampledir
    fdirs = os.listdir(fastq_path)
    fastq_path = rawdir+sampledir+'/'+''.join(fdirs)+'/'
    #print fastq_path
    image = imagedir+sampledir+'.tif'
    #print image
    resultdir2 = resultdir+sampledir+'/'
    pr.make_dir(resultdir2)
    job = "sp_cnt"
    slide = name[sampledir].split('-')[0] + '-' + name[sampledir].split('-')[1] 
    area = name[sampledir].split('-')[2].strip()
    samplename = sampledir
    #count_command = 'sbatch -J ' + job + ' -D ' + resultdir + ' --wrap="/home/jhy/spatial_transcriptomics/220905/spaceranger/spaceranger-2.0.0/spaceranger count --id=' + sampledir + ' --fastqs=' + fastq_path + ' --sample='+ samplename + ' --transcriptome=' + refdir + '"'
    #count_command = 'sbatch -J ' + job + ' -D ' + resultdir2 + ' --wrap="spaceranger count --id=' + sampledir + ' --fastqs=' + fastq_path + ' --sample='+ samplename + ' --transcriptome=' + refdir + ' --image=' + image  + '"'
    count_command = 'sbatch -J ' + job + ' -D ' + resultdir2 + ' --wrap="/home/jhy/spatial_transcriptomics/220905/spaceranger/spaceranger-2.0.0/spaceranger count --id=' + sampledir + ' --fastqs=' + fastq_path + ' --sample='+ samplename + ' --transcriptome=' + refdir + ' --image=' + image + ' --slide=' + slide + ' --area=' + area + '"'
    #print count_command
    pr.inShell(count_command)
    while (commands.getoutput('squeue -n' + job).count(job)>5): time.sleep(300)
    time.sleep(50)


