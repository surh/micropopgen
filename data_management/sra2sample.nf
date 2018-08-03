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
map = files(params.map)
sample_col = params.sample_col - 1
run_col  = params.run_col - 1

process process_run_list{
  cpus params.njobs

  input:
  file map from map

  output:
  set runs_per_sample into runs_per_sample

  exec:
  // Read map file line by line
  count = 0
  // map.eachLine{ str, i ->
  //   // Skip header
  //   if (count == 0 && params.header == true){
  //     count++
  //     return
  //   }
  //
  //   // Extract sample and run IDs
  //   def (sample, run) = str.split("\t")[sample_col, run_col]
  //
  //   // println "$sample\t$run"
  //   println tuple(sample, run)
  //   return tuple(sample, run)
  // }
  // .groupTuple()
  // .set{runs_per_sample}

  reader = map.newReader()
  while(str = reader.readLine()){
    // Skipe header
    if (count == 0 && params.header == true){
      count++
      continue
    }
    // Extract sample and run IDs
    (sample, run) = str.split("\t")[sample_col, run_col]

    // println "$sample\t$run"
    // println tuple(sample, run)
    Channel.from([tuple(sample, run)]).groupTuple().set{ runs_per_sample }
  }
}
