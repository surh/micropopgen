#!/usr/bin/env nextflow
// Copyright (C) 2021 Sur Herrera Paredes

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

// Nextflow pipeline that concatenates files according to a mapping
// file

// map must have two columns, tab-delimited, first column is path to
// files, second column name of new file
params.map = 'map.txt'
params.outdir = 'output'

FILES = Channel
    .fromPath(params.map)
    .splitCsv(header:false, sep:"\t")
    .map{ row -> tuple( file(row[0]), row[1] ) }
    .groupTuple(by:1)

process concatenate{
  publishDir params.outdir, mode: 'rellink'

  input:
  set file("infiles/"), newfile from FILES

  output:
  file newfile

  """
  cat infiles/* > $newfile
  """

}


// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns,owners'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '5h'
  memory = '1G'
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
