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

// Nextflow pipeline that takes output from metawas (GEMMA lmm results)
// and zipped genome feature files from midas db, and associates
// every SNP below a certain p-value threshold with a gene.

// Params
params.midas_db = ''
params.lmm_res = ''
params.genomes_file = ''
params.pval_thres = '1e-8'
params.outdir = 'output/'

// process params
genomes_file = file(params.genomes_file)

reader = genome_file.newReader()
GENOMES = []
while( line = reader.readLine() ) {
  GENOMES = GENOMES + [tuple(line,
    file("${params.lmm_res}/${line}_lmm.assoc.txt"),
    file("${params.midas_db}/rep_genomes/${line}/genome.features.gz"))]
}

process snps_to_genes{
  publisDir params.outdir, mode: 'copy', pattern: "*.closest"
  maxForks 10

  input:
  set genome, lmm_file, feat_file from GENOMES

  output:
  file 'snps.bed'
  file 'genome.features.bed'
  file "${genome}.closest"

  """
  # Convert snps to BED
  awk '(\$6 <= ${params.pval_thres}){print \$1 "\\t" \$2 "\\t" \$2}' \
    ${lmm_file} | sort -k1,1 -k2,2n | grep -v chr > snps.bed

  # Convert features to BED
  zcat ${feat_file} | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' | \
    grep -v scaffold_id | sort -k1,1 -k2,2n > genome.features.bed

  # Find closest
  closestBed -D a -a temp.snps.bed -b genome.features.bed > ${genome}.closest
  """
}



// 19637  awk '{print $1 $2 $2}'
// 19638  awk '{print $1 $2 $2}'  temp.snps
// 19639  awk '{print "$1\t$2 $2}'  temp.snps
// 19640  awk '{print "$1\t$2 $2"}'  temp.snps
// 19641  awk '{print $1 "\t" $2 "\t" $2}'  temp.snps
// 19642  awk '{print $1 "\t" $2 "\t" $2}'  temp.snps > temp.snps.bed
// 19644  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features
// 19645  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features > genome.features.bed
// 19653  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | sort -k1,1 -k2,2n
// 19654  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | grep -v scaffold_id | sort -k1,1 -k2,2n | \n\nls
// 19655  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | grep -v scaffold_id | sort -k1,1 -k2,2n\n\nls
// 19656  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | sort -k1,1 -k2,2n\n\nls
// 19660  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | sort -k1,1 -k2,2n\n\nls
// 19661  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | sort -k1,1 -k2,2n
// 19662  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | grep -v scaffold_id | sort -k1,1 -k2,2n
// 19663  awk '{print $2 "\t" $3 "\t" $4 "\t" $1}'  genome.features | grep -v scaffold_id | sort -k1,1 -k2,2n > genome.features.bed
// 19664  awk '{print $1 "\t" $2 "\t" $2}'  temp.snps | sort -k1,1 -k2,2n
// 19665  awk '{print $1 "\t" $2 "\t" $2}'  temp.snps | sort -k1,1 -k2,2n | head
// 19666  awk '{print $1 "\t" $2 "\t" $2}'  temp.snps | sort -k1,1 -k2,2n | grep -v chr
// 19667  awk '{print $1 "\t" $2 "\t" $2}'  temp.snps | sort -k1,1 -k2,2n | grep -v chr > temp.snps.bed
// 19683* awk '($6 < 1e-3)' output/*.txt | cut -f 1 | sort | uniq -c
