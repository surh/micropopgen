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

// Params
params.indir = ''
params.outdir = './output/'
params.dnds_mao = 'dnds.mao'
params.tajima_mao = 'tajima.mao'
params.njobs = 4

dnds_mao = file(params.dnds_mao)
tajima_mao = file(params.tajima_mao)

ALNS = Channel.fromFilePairs("${params.indir}/**.fasta",
  size: 1){
    file ->
      // gene_name = file.getName()
      // gene_name = gene_name.replaceAll(/.*fasta$/, '')
      // println gene_name
      file.name.replaceAll(/\.fasta$/, '')
    }
ALNS.into{ALNS_dnds; ALNS_tajima}
println "==============="

process mega_dnds{
  label 'mega'
  publishDir "${params.outdir}/dnds"
  maxForks params.njobs

  input:
  file dnds_mao
  set gene, file(aln) from ALNS_dnds

  output:
  file "${gene}.meg"

  """
  megacc -a $dnds_mao -d $aln -o ${gene}.dummysuffix
  """

}

process mega_tajima{
  label 'mega'
  publishDir "${params.outdir}/tajima"
  maxForks params.njobs

  input:
  file tajima_mao
  set gene, file(aln) from ALNS_tajima

  output:
  file "${gene}_summary.txt"

  """
  megacc -a $tajima_mao -d $aln -o ${gene}.dummysuffix
  """
}
