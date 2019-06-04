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

// Params
params.indir = ''
params.outdir = './output/'

// FEAT_FILES = Channel.fromPath("${params.indir}/*/genome.features.gz")
FEAT_FILES = Channel.fromFilePairs("${params.indir}/*/genome.features.gz",
  size: 1){
    file ->
      specname = file.getParent().getName()
      println specname
      // file.renameTo("${specname}.genome.features.gz")
      file.name.replaceAll(/genome.features.gz$/, "${specname}.closest.txt")}

println "================"

process closest_cds{
  publishDir params.outdir, saveAs: {"$spec"}

  input:
  set spec, file(features) from FEAT_FILES

  output:
  file "closest.txt"
  // 
  // exec:
  // println spec
  // println features

  // script:
  """
  $workflow.projectDir/find_closest_cds.sh $features
  """
}
