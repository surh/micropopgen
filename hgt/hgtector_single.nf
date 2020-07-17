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

params.indir = ''
params.close_tax = ''
params.genome_taxids = ''

genome_taxids = file(params.genome_taxids)

process read_genome_ids{
  label 'py3'

  input:
  file genome_taxids

  output:
  val x into OUT

  exec:

  myFile = file(genome_taxids)
  allLines  = myFile.readLines()
  for( line in allLines ) {
      println line
  }
  x = 1
}

OUT.subscribe{ println it }
