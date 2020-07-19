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
params.close_tax = ''
params.genome_taxids = ''

// Read tax ids from CSV
genome_taxids = file(params.genome_taxids)
TAXIDS = Channel
    .fromPath(genome_taxids)
    .splitCsv(header:true)
    .map{ row -> tuple(row.spec, row.tax_id) }

// Get list of input files
search_dir = file(params.search_dir)
SEARCHFILES = Channel.fromPath("$search_dir/*.tsv").
  .map{ search_file -> tuple(search_file.name.replaceAll(/\.tsv/, ""),
    file(search_file))}


process hgtector_analyse{
  label 'hgtector'

  input:
  tuple spec, file(search_file), taxid from SEARCHFILES.join(TAXIDS)

  exec:
  println "$spec\t$taxid\t$search_file"

}
