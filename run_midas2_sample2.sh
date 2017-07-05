#!/bin/bash
#PBS -N MIDAS2
#PBS -d ~/micropopgen/exp/2017/today4/
#PBS -o ~/micropopgen/exp/2017/today4/logs/midas2_sample2.log
#PBS -e logs/download.err
#PBS -l mem=2000mb
module load MIDAS
run_midas.py species midas2/sample1 -1 /godot/hmp/WGS/temp/test/SRS103943.tar.bz2
