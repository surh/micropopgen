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
params.genome_metadata = "metadata.txt"
params.min_size = 5

indir = file(params.indir)
snv_dir = file(params.snv_dir)
genome_metadata = file(params.genome_metadata)

// Species names that have SNVs in the catalogue
SPECSWITHSNVS = Channel.fromPath("$snv_dir/*_snvs.tsv")
  .map{snv_file -> tuple(snv_file.name.replaceAll(/_snvs\.tsv/, ''),
    file(snv_file))}
  // .into{SPECSWITHSNVS1; SPECSWITHSNVS2; SPECSWITHSNVS3; SPECSWITHSNVS4}

// Species that have a genome dir in the UHGG catalogue.
SPECSWITHGENOME = Channel.fromPath("$indir/*/*", type: 'dir')
  .map{specdir -> tuple(specdir.name, file(specdir))}
  // .into{SPECSWITHGENOME1; SPECSWITHGENOME2; SPECSWITHGENOME3; SPECSWITHGENOME4}

SPECSWITHSNVS.join(SPECSWITHGENOME)
  .into{INSPECS1; INSPECS2; INSPECS3; INSPECS4}

// Species to convert from tsv to vcf. The fasta file is
// needed to annotate contig sizes in the VCF header.
// UHGG2VCF = Channel.fromPath("$indir/*/*", type: 'dir')
//   .map{specdir -> tuple(specdir.name, file(specdir))}
//   .map{spec, specdir ->
//     tuple(spec,
//       file("$specdir/genome/${spec}.fna"),
//       file("$snv_dir/${spec}_snvs.tsv"))}
UHGG2VCF = INSPECS1
  .map{spec, snv_file, genome_dir ->
    tuple(spec, file("$genome_dir/genome/${spec}.fna"), file(snv_file))}
// UHGG2VCF.subscribe{ println it }

// FNA FILES to be split. Only the ones that have SNVs
// UHGGFNA = Channel.fromPath("$indir/*/*", type: 'dir')
//   .map{specdir -> tuple(specdir.name, file(specdir))}
//   .map{spec, specdir ->
//     tuple(spec,
//       file("$specdir/genome/${spec}.fna"))}
UHGGFNA = INSPECS2
  .map{spec, snv_file, genome_dir ->
    tuple(spec, file("$genome_dir/genome/${spec}.fna"))}

// GFF files for genomes that have SNVs. It is used to
// determine synonymous and non synonymous
// UHGGGFF = Channel.fromPath("$indir/*/*", type: 'dir')
//   .map{specdir -> tuple(specdir.name, file(specdir))}
//   .map{spec, specdir ->
//     tuple(spec,
//       file("$specdir/genome/${spec}.gff"))}
UHGGGFF = INSPECS3
  .map{spec, snv_file, genome_dir ->
    tuple(spec, file("$genome_dir/genome/${spec}.gff"))}

// SNV and GFF files to map SNVs to features. Only for genomes
// with SNVs
// SNVGFF = Channel.fromPath("$indir/*/*", type: 'dir')
//   .map{specdir -> tuple(specdir.name, file(specdir))}
//   .map{spec, specdir ->
//     tuple(spec,
//       file("$snv_dir/${spec}_snvs.tsv"),
//       file("$specdir/genome/${spec}.gff"))}
SNVGFF = INSPECS4
  .map{spec, snv_file, genome_dir ->
    tuple(spec,
      file(snv_file),
      file("$genome_dir/genome/${spec}.gff"))}

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

TABIXED.into{TABIXED1; TABIXED2}

process split_fnas{
  label 'py3'
  tag "$spec"
  // publishDir "$params.outdir/split_fnas", mode: 'rellink'

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

TABIXED_FNAS = TABIXED1.cross(SPLITFNAS1
  .map{spec, ctg_file -> tuple(spec,
    ctg_file.name.replaceAll(/\.fasta/, ""))})
  .map{vec1, vec2 -> tuple(vec2[0], vec2[1], file(vec1[1]), file(vec1[2]))}

process split_vcfs{
  label 'htslib'
  tag "${spec}.${ctg}"
  // publishDir "$params.outdir/ctg_vcfs", mode: 'rellink'

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
  // publishDir "$params.outdir/snv_effect", mode: 'rellink'

  input:
  tuple spec, ctg, file(vcf), file(gff), file(fasta) from SNVEFFIN

  output:
  tuple spec, file("${ctg}.tsv") optional true into SNVEFFS

  """
  mkdir vcf gff
  cd vcf
  ln -s ../$vcf ./
  cd ../
  awk '(\$1 == "$ctg")' $gff  | \
    sed "s/[']//g" > gff/${ctg}.gff

  Rscript ${workflow.projectDir}/snv_effect.r \
    vcf/ \
    gff/ \
    $fasta \
    ${ctg}.tsv
  """
}

process cat_snvs{
  label 'r'
  tag "$spec"
  publishDir "$params.outdir/snvs", mode: 'rellink'

  input:
  tuple spec, file("snveffs*.tsv") from SNVEFFS.groupTuple()

  output:
  tuple spec, file("${spec}.tsv") into SPECSNVEFFS

  """
  Rscript ${workflow.projectDir}/cat_tables.r \
    ${spec}.tsv \
    *.tsv
  """
}

process snvs2feats{
  label 'py3'
  tag "$spec"
  publishDir "$params.outdir/snv_feats", mode: 'rellink'

  input:
  tuple spec, file(snvs), file(gff) from SNVGFF

  output:
  tuple spec, file("${spec}.tsv") into SNV2FEATS

  """
  cut -f 1,2 $snvs | \
    grep -vP '^Contig' | \
    awk '{print \$1 "\t" \$2 - 1 "\t" \$2}' > snvs.bed

  cat $gff | \
    cut -f 1,3,4,5,9 | \
    sed 's/;/\t/' | \
    cut -f 1-5 | \
    sed 's/ID=//' | \
    awk '{print \$1 "\t" \$3 - 1 "\t" \$4 "\t" \$2 ";" \$5}' > genes.bed

  bedtools intersect \
    -wb \
    -a snvs.bed \
    -b genes.bed | \
    sort -k1V -k3n | \
    cut -f 1,3,7 | \
    sed 's/;/\t/' > ${spec}.tsv
  """
}

MKIN = TABIXED2
  .map{spec, vcf, tbi -> tuple(spec, file(vcf))}
  .join(SNV2FEATS)
  .join(SPECSNVEFFS)

process mktest{
  label 'r'
  label 'bigmem'
  label 'long'
  tag "$spec"
  publishDir "$params.outdir/mktest", mode: 'rellink'

  input:
  tuple spec, file(vcf),
    file("${spec}_snvfeats.tsv"),
    file("${spec}_snveffs.tsv") from MKIN
  file genome_metadata
  val min_size from params.min_size

  output:
  tuple spec, file("${spec}_mktest.tsv") optional true into MKTEST

  """
  Rscript ${workflow.projectDir}/mktest.r \
    ${spec}_snveffs.tsv \
    $vcf \
    ${spec}_snvfeats.tsv \
    $genome_metadata \
    ${spec}_mktest.tsv \
    $min_size
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 100
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
    module = 'R/3.6.1'
  }
  withLabel: 'long'{
    time = '24h'
  }
  withLabel: 'bigmem'{
    memory = '5G'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
