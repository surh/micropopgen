#!/usr/bin/env nextflow
// Copyright (C) 2018 Sur Herrera Paredes


params.genomes = 'genome_list.txt'
params.outdir = 'CDS'

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
  module 'eggnog'
  publishDir params.outdir
  time 24.h
  memory 2.GB

  input:
  val genome from GENOMES
  file faa from CDS

  output:
  file "${genome}.emapper.annotations" 

  exec:
  println genome

  script:
  """
  emapper.py --database bact --data_dir /opt/pkgs/eggnog/1.0.3/data/ --output_dir ./ -i ${faa} --cpu 1 --output ${genome}

  """

}


