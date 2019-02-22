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

// Nextflow pipeline that runs a directory of genomes through checkm.
// Combines results in single table

// Parameters:
// --indir
// Input directory with GENOME files. Assumes genomes are organized in
// an an arbitrary number of subdirectories. And each genome is in a
// second level subdirectory, where there is an '.fna' file with the
// same name as the genome directory. For example, the following structure
// would be acceptable:
// genomes/group1/genome1/genome1.fna
// genomes/group1/genome2/genome2.fna
// genomes/group2/genome3/genome3.fna
// genomes/group2/genome4/genome4.fna

// --outdir
// Name of output directory

// --conda_checkm
// Path to conda environment where checkm is installed. Default for fraserv

// --batch_size
// Number of genomes to analyze by batch

// --threads, --memory, --time, --queue, --max_forks
// Parameters for checkm jobs.


// Params
params.indir = 'genomes'
// params.outfile = 'checkm_results.txt'
params.outdir = 'output/'
params.contamination = 2
params.completeness = 98

params.conda_checkm = '/opt/modules/pkgs/anaconda/3.6/envs/python2'
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
  conda params.conda_checkm
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

  # tail -n +2 checkm_results.txt > checkm_results_noheader.txt
  """
}

process collect_results{
  cpus 1
  memory '1GB'
  time '00:30:00'
  module 'fraserconda'
  queue params.queue
  publishDir params.outdir

  input:
  file "*.txt" from checkm_results.collect()

  output:
  file "checkm_results.txt" into CHECKM

  """
  ${workflow.projectDir}/../sutilspy/bin/cat_tables.py \
    *.txt \
    --outfile checkm_results.txt
  """
}

process filter_checkm{
  cpus 1
  memory '2GB'
  time '00:30:00'
  queue params.queue
  publishDir "${params.outdir}/checkm}"

  input:
  file "checkm_results.txt" from CHECKM

  output:
  file "output/all_completeness_histogram.svg"
  file "output/all_contamination_histogram.svg"
  file "output/all_heterogeneity_histogram.svg"
  file "output/all_completeness_vs_contamination.svg"
  file "output/chosen_completeness_histogram.svg"
  file "output/chosen_contamination_histogram.svg"
  file "output/chosen_heterogeneity_histogram.svg"
  file "output/chosen_completeness_vs_contamination.svg"
  file "output/chosen_checkm_results.txt" into CHOSEN

  """
  ${workflow.projectDir}/process_checkm_results.r \
    checkm_results.txt \
    --contamination ${params.contamination} \
    --completeness ${params.completeness}
  """
}
