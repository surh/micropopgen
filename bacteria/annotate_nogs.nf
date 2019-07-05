#!/usr/bin/env nextflow
// Copyright (C) 2019 Sur Herrera Paredes

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

// Params
params.indir = ''
params.db = ''
params.outdir = 'output'

FILES = Channel.fromPath("${params.indir}/*")

process annotatate_nogs{
  label 'py3'
  publishDir params.outdir, saveAs: {"$input"}, mode: 'rellink'

  input:
  file input from FILES

  output:
  file 'eggnog_table.txt'

  """
  ${workflow.projectDir}/annotate_nogs.py \
    --input $input \
    --db ${params.db}
  """
}

// example nextflow.log
/*
process{
  maxForks = 30
  withLabel: py3{
    module = 'fraserconda'
  }
}
exectuor{
  name = 'slurm'
  queueSize = 100
}
*/
