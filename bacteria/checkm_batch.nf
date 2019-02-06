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
params.outfile = 'checkm_results.txt'
params.batch_size = 200
params.threads = 8
params.memory = '40GB'
params.time = '2:00:00'
params.queue = 'owners,hbfraser,hns'
params.max_forks = 200

// Process params
indir = file(params.indir)

process create_batch_map{
  cpus 1
  memory '1GB'
  time '1:00:00'
  errorStrategy 'retry'
  maxRetries 2
  queue params.queue
  module 'fraserconda'

  input:
  file indir

  output:
  file('batch_map.txt') into create_batch
  file('checkm_batches/batch_*') into checkm_dirs

  """
  ${workflow.projectDir}/create_batches.py --indir ${params.indir} \
    --outdir checkm_batches \
    --outfile batch_map.txt \
    --batch_size ${params.batch_size}
  """
}

process run_checkm{
  cpus params.threads
  memory params.memory
  time params.time
  errorStrategy 'retry'
  maxRetries 2
  maxForks params.max_forks
  module 'prodigal:hmmer:pplacer:fraserconda'
  conda 'python2'
  queue params.queue

  input:
  file checkm_dir from checkm_dirs.flatten()

  output:
  file "checkm_results.txt" into checkm_results
  // file "checkm_results_noheader.txt" into checkm_results_noheader

  """
  checkm lineage_wf \
    -t ${params.threads} \
    -f checkm_results.txt \
    --tab_table \
    ${checkm_dir} \
    results

  // tail -n +2 checkm_results.txt > checkm_results_noheader.txt
  """
}

process collect_results{
  cpus 1
  memory '1GB'
  time '00:30:00'
  module 'fraserconda'
  queue params.queue
  publishDir './'

  input:
  file '*.txt' from checkm_results.collect()

  output:
  file "${params.outfile}"

  """
  ${workflow.projectDir}/../sutilspy/bin/cat_tables.py \
    *.txt \
    --outfile ${params.outfile}
  """
}

// checkm_results_noheader
//   .collectFile()
//   .println{ it.text }
