#!/usr/bin/env nextflow
// Copyright (C) 2018-2019 Sur Herrera Paredes


params.genomes = 'genomes.txt'
params.outdir = 'annots'
params.njobs = 20

genomes = file(params.genomes)
reader = genomes.newReader()
GENOMES = []
CDS = []
while(str = reader.readLine()){
  GENOMES = GENOMES + [str]
  cds = file("/godot/users/sur/data/genomes/midas_db_v1.2/FAA/${str}.CDS.faa")
  CDS = CDS + [cds]
}

process eggnog{
  label 'eggnog'
  publishDir params.outdir
  maxForks params.njobs

  input:
  val genome from GENOMES
  file faa from CDS

  output:
  file "${genome}.emapper.annotations" 

  exec:
  println genome

  script:
  """
  emapper.py --database bact \
    --data_dir /opt/pkgs/eggnog/1.0.3/data/ \
    --output_dir ./ \
    -i ${faa} \
    --cpu 1 \
    --output ${genome}
  """
}


// Example nextflow.config
/*
process {
  executor = 'slurm'
  withLabel: eggnog{
    module = 'eggnog'
    time = '24h'
    memory = '2G'
  }
}
*/
