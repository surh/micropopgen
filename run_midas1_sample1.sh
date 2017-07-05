#!/bin/bash
#PBS -N MIDAS1.1
#PBS -d /home/sur/micropopgen/exp/2017/today4/
#PBS -o logs/midas1_sample1.log
#PBS -e logs/midas1_sample1.err
#PBS -l mem=2000mb
module load MIDAS
pwd
run_midas.py species midas1/sample1 -1 /godot/hmp/WGS/temp/temp_fastq/SRS011061_R1.fastq.gz  -2 /godot/hmp/WGS/temp/temp_fastq/SRS011061_R2.fastq.gz
