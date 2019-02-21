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

// Nextflow pipeline that splits patric dowloads to avoid being disconnected.

params.genomes = 'genomes.txt'
params.chunk_size = '300'
params.id_col = 1
params.group_col = 5
params.name_col = 0
params.failed = 'failed.txt'
params.outdir = 'patric'
params.max_forks = 2


genomes = file(params.genomes)

process split_genomes{
  module 'fraserconda'

  input:
  file genomes

  output:
  file 'tab*.txt' into CHUNKS

  """
  ${workflow.projectDir}/../sutilspy/bin/split_tables.py \
    ${genomes} \
    --nlines ${params.chunk_size}
  """
}

process download{
  memory '2GB'
  time '5:00:00'
  module 'fraserconda'
  publishDir params.outdir, pattern: 'patric/*'
  maxForks params.max_forks

  input:
  file 'genomes_chunk.txt' from CHUNKS.flatten()

  output:
  file 'failed.txt' optional true into FAILED
  file 'patric/*'


  """
  ${workflow.projectDir}/patric_download_genomes.py \
    --genomes genomes_chunk.txt \
    --outdir patric \
    --check \
    --id_col ${params.id_col} \
    --group_col ${params.group_col} \
    --name_col ${params.name_col} \
    --header \
    --failed failed.txt
  """
  // ~/micropopgen/src/micropopgen/patric/patric_download_genomes.py
  // --genomes selected_strains3.txt --outdir patric --check
  // --id_col 1 --group_col 5 --name_col 0
  // --header --failed failed.txt &> dowload_patric3.log
}

process collect_failed{
  cpus 1
  memory '1GB'
  time '00:30:00'
  module 'fraserconda'
  publishDir './'

  input:
  file failed_genomes from FAILED.collect()

  output:
  file "${params.outfile}"

  """
  ${workflow.projectDir}/../sutilspy/bin/cat_tables.py \
    $failed_downloads \
    --outfile ${params.failed}
  """
}
