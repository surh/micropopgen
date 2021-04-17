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

// Params
params.snv_dir = ''
params.uhgg_catalogue = ''
params.map = ''
params.min_genomes = 5
params.only_isolates = false
params.outdir  = 'output'

// Process params,
snv_dir = file(params.snv_dir)
uhgg_catalogue = file(params.uhgg_catalogue)
map = file(params.map)
only_isolates = ''
if (params.only_isolates){
  only_isolates = '--only_isolates'
}

SNVS = Channel.fromPath("$snv_dir/*_snvs.tsv")
  .map{ snvfile -> tuple(snvfile.name.replaceAll(/_snvs\.tsv$/, ''),
    file(snvfile)) }

GFFS = Channel.fromPath("$uhgg_catalogue/*/*/genome/*.gff")
  .map{ gfffile -> tuple(gfffile.name.replaceAll(/\.gff$/, ''),
    file(gfffile)) }


process fst{
  tag "$spec"
  label 'r'
  publishDir "$params.outdir/n_genomes", mode: 'rellink',
    pattern: 'output/n_genomes.tsv', saveAs: {"${spec}.tsv"}
  publishDir "$params.outdir/site_fst", mode: 'rellink',
    pattern: 'output/site_fst.tsv', saveAs: {"${spec}.tsv"}
  publishDir "$params.outdir/feature_fst", mode: 'rellink',
    pattern: 'output/feature_fst.tsv', saveAs: {"${spec}.tsv"}

  input:
  tuple spec, file(snvs), file(gff) from SNVS.join(GFFS)
  val min_genomes from params.min_genomes
  val only_isolates from only_isolates
  file map from map

  output:
  file "output/n_genomes.tsv"
  file "output/site_fst.tsv"
  file "output/feature_fst.tsv" optional true

  """
  Rscript $workflow.projectDir/all_genes_fst.r \
    $snvs \
    $map \
    --gff $gff \
    --pop_col Country \
    --pops "United States" "China" \
    $only_isolates \
    --outdir output \
    --min_genomes $min_genomes
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '72h'
  memory = '10G'
  withLabel: 'r'{
    module = 'R/3.6.1'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
