#!/bin/bash
#PBS -N MIDAS2.2
#PBS -d /home/sur/micropopgen/exp/2017/today4/
#PBS -o logs/midas2_sample1.log
#PBS -e logs/midas2_sample1.err
#PBS -l mem=2000mb
module load MIDAS
pwd
run_midas.py species midas2/sample1 -1 /godot/hmp/WGS/temp/test/SRS013216/SRS013216.fastq
