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


indir = file(params.indir)
snv_dir = file(params.snv_dir)

UHGG2VCF = Channel.fromPath("$indir/*/*", type: 'dir')
  .map{specdir -> tuple(specdir.name, file(specdir))}
  .map(spec, specdir <- tuple(spec,
    file("$specdir/genome/${spec}.fna"),
    file("$snv_dir/${spec}_snvs.tsv")))

UHGG2VCF.subscribe{ println it }
