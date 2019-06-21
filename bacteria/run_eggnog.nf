#!/usr/bin/env nextflow
// Copyright (C) 2018-2019 Sur Herrera Paredes

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

// Nextflow pipeline that runs a directory of genomes through checkm.
// Combines results in single table


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
