#!/usr/bin/env nextflow
// Copyright (C) 2017 Sur Herrera Paredes

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Nextflow pipeline that goes from genome files as downloaded
// from PATRIC, formats the data for roary, and runs roary to
// obtain a set of orthologous groups.


// Main parameters
params.indir = 'runs/'
params.fastq_dir = 'fastq/'
params.outdir = 'samples/'
params.map = 'map.txt'

// Other Parameters
params.sample_col = 1
params.run_col = 2
params.keep_intermediate = false
params.header = true
params.failed = 'failed.txt'
params.njobs = 10
params.failed_dir = 'failed/'
params.submissions_dir = 'submissions/'
params.logdir = 'logs/'

// Process params
map = file(params.map)
sample_col = params.sample_col - 1
run_col  = params.run_col - 1

// Read mapping file
count = 0
reader = map.newReader()
run_sample_table = []
while(str = reader.readLine()){
  // Skip header
  if (count == 0 && params.header == true){
    count++
    continue
  }
  // Extract sample and run IDs
  (sample, run) = str.split("\t")[sample_col, run_col]
  run_sample_table = run_sample_table + [tuple(sample, run)]
}
// Convert to channel that maps sample to its runs
Channel.from(run_sample_table)
  .map{sample, run ->
    return tuple(sample, run)}
  .groupTuple()
  .set{runs_groups}
// Get channel with runs
Channel.from(run_sample_table).map{sample, run -> return tuple(run, file("${params.indir}/${run}.sra"))}.set{runs}


// Convert all sra files to fastq
process fastqdump{
  cpus 1
  maxForks params.njobs

  input:
  set run, run_file from runs

  output:
  set file("${params.fastq_dir}/${run}_1.fastq.bz2"),
    file("${params.fastq_dir}/${run}_1.fastq.bz2") into fastq_files

  """
  # cat ${run_file}.sra
  fastq-dump -I -O ${params.fastq_dir} --split-files --bzip2 ${run_file}
  """
}

// process sra2fastq{
//   cpus 1
//   maxForks params.njobs
//
//   input:
//   set sample, runs from runs_groups
//
//   exec:
//   println sample + runs
//
//   // Within each group of runs corresponding to the same sample,
//   // use fastq-dump to convert them
//   process fastqdump{
//     cpus 1
//     maxForks params.njobs
//
//
//     input:
//     val run from runs
//     val sample from sample
//     file "${params.indir}/${run}.sra" from runs
//
//     output:
//     set sample, run, file("${params.fastq_dir}/${run}_1.fastq.bz2"), file("${params.fastq_dir}/${run}_1.fastq.bz2") into fastq_files
//
//     """
//     # cat ${params.indir}/${run}.sra
//     fastq-dump -I -O ${params.fastq_dir} --split-files --bzip2 ${params.indir}/${run}.sra
//     """
//   }
//
// }

//
// println runs_per_sample
// runs_per_sample.groupTuple().set{ run_groups }
//
//
// process fastq_dump_runs{
//   cpus params.njobs
//
//   input:
//   set sample, runs from run_groups
//
//   exec:
//   println sample
// }
//
// process validate_run{
//   cpus params.njobs
//
//   input:
//   file "params.indir/$run\.sra" from runs
//
//   """
//   vdb-validate params.indir/$run\.sra
//   """
// }
