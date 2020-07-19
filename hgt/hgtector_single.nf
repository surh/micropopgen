#!/usr/bin/env nextflow
// Copyright (C) 2020 Sur Herrera Paredes

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

// Pipeline to run hgtector on a number of genomes individually

params.search_dir = ''
params.close_tax = -1
params.genome_taxids = ''
params.taxdump_dir = ''
params.outdir = "output/"

// Read tax ids from CSV
genome_taxids = file(params.genome_taxids)
TAXIDS = Channel
  .fromPath(genome_taxids)
  .splitCsv(header:true, sep:"\t")
  .map{ row -> tuple(row.spec,
    row.tax_id,
    (params.close_tax == -1) ? row.tax_id : params.close_tax) }.
    subscribe{ println it }

// Get list of input files
// search_dir = file(params.search_dir)
// SEARCHFILES = Channel.fromPath("$search_dir/*.tsv")
//   .map{ search_file -> tuple(search_file.name.replaceAll(/\.tsv/, ""),
//     file(search_file))}
//
// // Stage db dir
// taxdump_dir = file(params.taxdump_dir)
//
// process hgtector_analyse{
//   label 'hgtector'
//   tag "$spec"
//   publishDir params.outdir, mode: 'rellink'
//
//   input:
//   tuple spec, file(search_file), taxid, close_tax from SEARCHFILES.join(TAXIDS)
//   file taxdump_dir
//
//   output:
//   tuple spec, file("$spec/")
//
//   """
//   hgtector analyze \
//     --input $search_file \
//     --output $spec \
//     --taxdump $taxdump_dir \
//     --self-tax $taxid \
//     --close-tax $close_tax
//   """
//
// }


// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '48h'
  memory = '1G'
  withLabel: 'hgtector'{
    module = 'anaconda'
    conda = "/opt/modules/pkgs/anaconda/4.8/envs/hgtector"
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
