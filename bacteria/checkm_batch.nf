#!/usr/bin/env nextflow
// Copyright (C) 2018 Sur Herrera Paredes

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

// Nextflow pipeline that runs a directory of genomes through checkm

params.indir = 'genomes'
params.batch_size = 200
params.threads = 8
params.memory = '20GB'


process create_batch_map{
  output:
  file('batch_map.txt')

  """
  create_batch_map.py --indir ${params.indir} \
    --outfile batch_map.txt
    --batch_size ${params.batch_size}
  """
}
