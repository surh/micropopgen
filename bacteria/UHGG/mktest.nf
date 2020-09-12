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

UHGGGFF = Channel.fromPath("$indir/*/*", type: 'dir')
  .map{specdir -> tuple(specdir.name, file(specdir))}
  .map{spec, specdir ->
    tuple(spec,
      file("$specdir/genome/${spec}.gff"))}

process snvs2vcf{
  label 'py3'
  tag "$spec"

  input:
  tuple spec, file(genome_fna), file(snvs) from UHGG2VCF

  output:
  tuple spec, file("${spec}.vcf") into UNSORTEDVCF

  """
  ${workflow.projectDir}/uhgg_snv2vcf.py \
    --input $snvs \
    --genome_fasta $genome_fna \
    --output ${spec}.vcf \
    --include_genomes
  """
}

process tabix_vcf{
  label 'htslib'
  tag "$spec"
  publishDir "$params.outdir/tabix", mode: 'rellink'

  input:
  tuple spec, file(vcf) from UNSORTEDVCF

  output:
  tuple spec, file("${spec}.vcf.gz"), file("${spec}.vcf.gz.tbi") into TABIXED

  """
  (grep ^"#" $vcf; grep -v ^"#" $vcf | sort -k1,1V -k2,2n) | \
    bgzip > ${spec}.vcf.gz

  tabix -p vcf ${spec}.vcf.gz
  """
}

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
  .into{SPLITFNAS1; SPLITFNAS2}

// SPLITFNAS1
//   .map{spec, ctg_file -> tuple(spec,
//     ctg_file.name.replaceAll(/\.fasta/, ""),
//     file(ctg_file))}
//

TABIXED_FNAS = TABIXED.cross(SPLITFNAS1
  .map{spec, ctg_file -> tuple(spec,
    ctg_file.name.replaceAll(/\.fasta/, ""))})
  .map{vec1, vec2 -> tuple(vec2[0], vec2[1], file(vec1[1]), file(vec1[2]))}

process split_vcfs{
  label 'htslib'
  tag "${spec}.${ctg}"
  publishDir "$params.outdir/ctg_vcfs", mode: 'rellink'

  input:
  tuple spec, ctg, file(vcf), file(tbi) from TABIXED_FNAS

  output:
  tuple spec, ctg, file("${ctg}.vcf") into CTGVCF

  """
  zcat $vcf | grep -P '^#' > header.txt
  tabix $vcf $ctg > snvs.vcf

  cat header.txt snvs.vcf > ${ctg}.vcf
  """
}

// Prepare input for PopGenome script
VCF_GFFS = UHGGGFF
  .cross(CTGVCF)
  .map{vec1, vec2 -> tuple(vec2[0], vec2[1], file(vec2[2]), file(vec1[1]))}
  // .subscribe{println it}
// [MGYG-HGUT-00001, GUT_GENOME000001_90, /cashew/users/sur/exp/fraserv/2020/today/work/ca/dfeb29e0cab79ef333cccee8001fa5/GUT_GENOME000001_90.vcf, /cashew/users/sur/exp/fraserv/2020/today/genomes/MGYG-HGUT-000/MGYG-HGUT-00001/genome/MGYG-HGUT-00001.gff]
SPLITFNAFILES = SPLITFNAS2
  .map{spec, ctg_file -> tuple(spec,
    ctg_file.name.replaceAll(/\.fasta/, ""),
    file(ctg_file))}
  // .subscribe{println it}
// [MGYG-HGUT-00001, GUT_GENOME000001_96, /cashew/users/sur/exp/fraserv/2020/today/work/64/485951511c4158aadba8bb1c0f12ba/MGYG-HGUT-00001/GUT_GENOME000001_96.fasta]
SNVEFFIN = VCF_GFFS.join(SPLITFNAFILES, by: [0,1])

process snv_effect{
  label 'r'
  tag "${spec}.${ctg}"
  publishDir "$params.outdir/snv_effect", mode: 'rellink'

  input:
  tuple spec, ctg, file(vcf), file(gff), file(fasta) from SNVEFFIN

  output:
  tuple spec, ctg, file("${ctg}.tsv") optional true into SNVEFFS

  """
  mkdir vcf gff
  ln -s $vcf vcf/
  awk '(\$1 == "$ctg")' $gff > gff/${ctg}.gff

  Rscript ${workflow.projectDir}/snv_effect.r \
    vcf/ \
    gff/ \
    $fasta \
    ${ctg}.tsv
  """
}


// awk '($1 == "GUT_GENOME000004_1")' ../../gff/snvs_genomes.gff > snvs.gff


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
  withLabel: 'r'{
    module = 'R/4.0.2'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
