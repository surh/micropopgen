#!/bin/bash
#PBS -N MIDAS2
#PBS -o logs/midas2_sample2.log
#PBS -e logs/midas2_sample2.err
#PBS -l mem=2000mb
module load MIDAS
run_midas.py species midas2/sample1 -1 /godot/hmp/WGS/temp/test/SRS103943.tar.bz2
