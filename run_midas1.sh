#!/bin/bash
#PBS -N MIDAS1
#PBS -d ~/micropopgen/experiments/2017/today4/
#PBS -o ~/micropopgen/experiments/2017/today4/logs/download.log
#PBS -e logs/download.err
#PBS -l mem=2000mb
module load MIDAS
run_midas.py species midas1/ -1 /godot/hmp/WGS/temp/temp_fastq/SRS011061_R1.fastq.gz  -2 /godot hmp/WGS/temp/temp_fastq/SRS011061_R2.fastq.gz