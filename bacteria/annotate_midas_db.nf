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

// Take a directory with genomes in midas_db format. Extract all CDS, in fna
// translate to faa, and run eggnog on all of them.

// Params
params.midas_dir = ''
params.njobs = 20
params.outdir = 'output'
params.eggnog_threads = 4
params.eggnog_db = 'bact'

GENOMEDIRS = Channel.fromPath("${params.midas_dir}/*", type: 'dir', maxDepth: 0)

process extract_fna {
  label 'bedtools'
  maxForks params.njobs
  publishDir "${params.outdir}/FNA", mode: 'rellink'

  input:
  file spec from GENOMEDIRS

  output:
  file "${spec.getName()}.CDS.fna" into FNAS
  val "$spec" into SPECS_TRANSLATE, SPECS_EGGNOG

  """
  # Convert CDS features to BED
  zcat $spec/genome.features.gz | \
    awk '(\$6 == "CDS")' | \
    awk '{print \$2 "\t" \$3 - 1 "\t" \$4 "\t" \$1 "\t.\t" \$5}' | \
    sort -k1,1 -k2,2n > genome.features.bed

  # Uncompress genome
  zcat $spec/genome.fna.gz > genome.fna

  # Create fasta file
  bedtools getfasta \
    -fi genome.fna \
    -bed genome.features.bed \
    -s \
    -name | \
    sed 's/([\\+\\-])//' > ${spec.getName()}.CDS.fna

  # Clean
  rm genome.features.bed genome.fna
  """
}

process translate{
  label 'py3'
  maxForks params.njobs
  publishDir "${params.outdir}/FAA", mode: 'rellink'

  input:
  file fna_file from FNAS
  val spec from SPECS_TRANSLATE

  output:
  file "${spec}.CDS.faa" into FAAS

  """
  python ${workflow.projectDir}/translate.py \
    --infile $fna_file \
    --remove_stops \
    --outfile ${spec}.CDS.faa
  """
}

process eggnog{
  label 'eggnog'
  maxForks params.njobs
  publishDir "params.outdir/eggnog", mode: 'rellink'
  cpus params.eggnog_threads

  input:
  file faa_file from FAAS
  val spec from SPECS_EGGNOG

  output:
  file "${spec}.emapper.annotations" into ANNOTS

  exec:
  println genome

  script:
  """
  emapper.py \
    --database ${params.eggnog_db} \
    --data_dir $DATA_PATH \
    --output_dir ./ \
    -i $faa_file \
    --cpu ${params.eggnog_threads} \
    --output $spec
  """
}

// example nextflow.config
/*
process {
  executor = 'slurm'
  queue = 'hbfraser'
  label: 'bedtools' {
    module = 'anaconda'
    conda = 'fraserconda'
    cpu = 1
    time = '1h'
    memory = '1G'
  }
  label: 'py3' {
    module = 'anaconda'
    conda = 'fraserconda'
    cpu = 1
    time = '1h'
    memory = '1G'
  }
  label: 'eggnog'{
    module = 'eggnog'
    time = '24h'
    memory = '2G'
  }
}
*/
