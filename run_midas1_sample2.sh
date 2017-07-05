#!/bin/bash
#PBS -N MIDAS1
#PBS -d ~/micropopgen/exp/2017/today4/
#PBS -o ~/micropopgen/exp/2017/today4/logs/midas1_sample2.log
#PBS -e logs/download.err
#PBS -l mem=2000mb
module load MIDAS
run_midas.py species midas1/sample2 -1 /godot/hmp/WGS/temp/temp_fastq/SRS011098_R1.fastq.gz  -2 /godot/hmp/WGS/temp/temp_fastq/SRS011098_R2.fastq.gz
