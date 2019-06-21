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
params.indir = ''
params.njobs = 20
params.outdir = 'output'

indir = file(params.indir)
println indir
GENOMEDIRS = Channel.fromPath("${params.indir}/*", type: 'dir', maxDepth: 0)
process extract_fna {
  label 'bedtools'
  maxForks params.njobs
  publishDir "${params.outdir}/FNA", mode: 'rellink'

  input:
  file spec from GENOMEDIRS

  output:
  file "${spec.name}.CDS.fna" into FNAS

  """
  # Convert features to BED
  zcat $spec/genome.features.gz | \
    grep -v scaffold_id | \
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
    sed 's/([\\+\\-])//' > ${spec.name}.CDS.fna

  # Clean
  rm genome.features.bed genome.fna
  """
}
