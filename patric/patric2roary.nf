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

// Nextflow pipeline that goes from genome files as downloaded
// from PATRIC, formats the data for roary, and runs roary to
// obtain a set of orthologous groups.


// Parameters
params.indir = "patric/"
params.files = 'files.txt'
params.threads = 8
params.outdir = 'roary/'
params.njobs = 10

// Get list of files
files = file(params.files)



// PROCESSES
process preprocess_patric_gff{
  maxForks params.njobs

  input:
  file patric_gff from patric_gffs

  output:
  file 'roary.gff' into roary_gffs

  """
  ${workflow.projectDir}/patric2roary_gff.py \
    --infile ${patric_gff} \
    --outfile roary.gff
  """
}

// process preprocesss_patric_fna{
//   maxForks params.njobs
//
//   input:
//   file patric_fna from patric_fnas
//
//   output:
//   file 'roary.fna' into roary_fnas
//
//   """
//   ${params.bindir}/fasta_remove_desc.py \
//     --infile ${patric_fna} \
//     --outfile roary.fna
//   """
// }

process create_roary_input{
  maxForks params.njobs

  input:
  file gff from roary_gffs
  file fna from roary_fnas

  output:
  file 'roary_input.gff' into roary_inputs

  """
  cat ${gff} ${fna} > roary_input.gff
  """
}

process run_roary{
  cpus params.threads
  publishDir params.outdir, mode: 'move'

  input:
  file '*.gff' from roary_inputs.collect()

  output:
  file params.outdir

  """
  roary -p ${params.threads} \
    -f ${params.outdir} \
    -v \
    *.gff
  """
}
