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

// Params
params.indir = "genomes/"
params.snv_dir = "snv_catalogue/"
params.outdir = "output/"

indir = file(params.indir)
snv_dir = file(params.snv_dir)

UHGG2VCF = Channel.fromPath("$indir/*/*", type: 'dir')
  .map{specdir -> tuple(specdir.name, file(specdir))}
  .map{spec, specdir ->
    tuple(spec,
      file("$specdir/genome/${spec}.fna"),
      file("$snv_dir/${spec}_snvs.tsv"))}
// UHGG2VCF.subscribe{ println it }
UHGGFNA = Channel.fromPath("$indir/*/*", type: 'dir')
  .map{specdir -> tuple(specdir.name, file(specdir))}
  .map{spec, specdir ->
    tuple(spec,
      file("$specdir/genome/${spec}.fna"))}


// process snvs2vcf{
//   label 'py3'
//   tag "$spec"
//
//   input:
//   tuple spec, file(genome_fna), file(snvs) from UHGG2VCF
//
//   output:
//   tuple spec, file("${spec}.vcf") into UNSORTEDVCF
//
//   """
//   ${workflow.projectDir}/uhgg_snv2vcf.py \
//     --input $snvs \
//     --genome_fasta $genome_fna \
//     --output ${spec}.vcf \
//     --include_genomes
//   """
// }

// process tabix_vcf{
//   label 'htslib'
//   tag "$spec"
//   publishDir "$params.outdir/tabix", mode: 'rellink'
//
//   input:
//   tuple spec, file(vcf) from UNSORTEDVCF
//
//   output:
//   tuple spec, file("${spec}.vcf.gz"), file("${spec}.vcf.gz.tbi")
//
//   """
//   (grep ^"#" $vcf; grep -v ^"#" $vcf | sort -k1,1V -k2,2n) | \
//     bgzip > ${spec}.vcf.gz
//
//   tabix -p vcf ${spec}.vcf.gz
//   """
// }

process split_fnas{
  label 'py3'
  tag "$spec"
  publishDir "$params.outdir/split_fnas", mode: 'rellink'

  input:
  tuple spec, file(genome_fna) from UHGGFNA

  output:
  tuple spec, file("$spec/*.fasta") into SPLITFNAS

  """
  ${workflow.projectDir}/extract_contigs.py \
    --input $genome_fna \
    --outdir $spec
  """
}

SPLITFNAS
  .transpose()
  .subscribe{println it}
  // .transpose()
  // .view()


// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 40
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '5h'
  memory = '1G'
  withLabel: 'py3'{
    module = 'anaconda'
    conda = "/opt/modules/pkgs/anaconda/4.8/envs/fraserconda"
  }
  withLabel: 'htslib'{
    module = 'htslib'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
