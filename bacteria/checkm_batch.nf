#!/usr/bin/env nextflow
// Copyright (C) 2018 Sur Herrera Paredes

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

// Nextflow pipeline that runs a directory of genomes through checkm

params.indir = 'genomes'
params.batch_size = 200
params.threads = 8
params.memory = '20GB'


process create_batch_map{
  cpus 1
  memory '1GB'
  time 1:00:00
  errorStrategy 'retry'
  maxRetries 3

  output:
  file('batch_map.txt') into create_batch
  file('checkm_batches/batch_*') into checkm_dirs

  """
  create_batches.py --indir ${params.indir} \
    --outdir checkm_batches \
    --outfile batch_map.txt \
    --batch_size ${params.batch_size}
  """
}

process run_checkm{
  cpus params.threads
  memory '10GB'
  time '2:00:00'
  errorStrategy 'retry'
  maxRetries 3
  maxForks 200
  module 'prodigal:hmmer:pplacer:fraserconda'
  conda '/share/PI/hbfraser/modules/packages/anaconda3/5.1/envs/python2/'

  input:
  file checkm_dir from checkm_dirs

  """
  checkm lineage_wf \
    -t ${params.threads} \
    -f checkm_results.txt \
    --tab_table genomes \
    ${checkm_dir} \
    results
  """

}
