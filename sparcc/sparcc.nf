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

/*
34419  python ../SparCC.py fake_data.txt -i 5
34422  python ../SparCC.py fake_data.txt -i 5
34423* cd /opt/modules/pkgs/anaconda/3.6/envs/sparcc
34429  python ../SparCC.py fake_data.txt -i 5
34431  more cor_mat_SparCC.out

34432  python ../MakeBootstraps.py -h
34433  python ../MakeBootstraps.py fake_data.txt -n 5 -t permutation_#.txt -p pvals/

34438  python ../SparCC.py pvals/permutation_1.txt -i 5 --cor_file=pvals/perm_cor_1.txt
34439  python ../SparCC.py pvals/permutation_0.txt -i 5 --cor_file=pvals/perm_cor_0.txt
34440  python ../SparCC.py pvals/permutation_2.txt -i 5 --cor_file=pvals/perm_cor_2.txt
34441  python ../SparCC.py pvals/permutation_3.txt -i 5 --cor_file=pvals/perm_cor_3.txt
34442  python ../SparCC.py pvals/permutation_4.txt -i 5 --cor_file=pvals/perm_cor_4.txt
34445  python ../PseudoPvals.py cor_mat_SparCC.out pvals/permutation_#.txt
34446  python ../PseudoPvals.py cor_mat_SparCC.out pvals/permutation_#.txt 5
34447  python ../PseudoPvals.py cor_mat_SparCC.out pvals/permutation_#.txt 5 -h
34448  python ../PseudoPvals.py cor_mat_SparCC.out pvals/permutation_#.txt 5 -h -o pvals.txt
34449  python ../PseudoPvals.py cor_mat_SparCC.out pvals/permutation_#.txt 5 -o pvals.txt
34453  python ../PseudoPvals.py cor_mat_SparCC.out pvals/perm_cor_#.txt 5 -o pvals.txt
*/

// Params
params.input = ''
params.iter = 20
params.perms = 100
params.outdir = 'output'

input = file(params.input)

process sparcc_cor{
  label 'sparcc'
  publishDir "${params.outdir}/cor", mode: 'rellink'

  input:
  file input

  output:
  file 'cor_mat_SparCC.out' into COR
  file 'cov_mat_SparCC.out'

  """
  SparCC.py $input -i ${params.iter}
  """
}

process sparcc_bootstraps{
  label 'sparcc'
  publishDir "${params.outdir}", mode: 'rellink'

  input:
  file input

  output:
  file 'perms/permutation_*.txt' into PERMS

  """
  MakeBootraps.py $input -n ${params.perms} -t permutation_#.txt -p perms/
  """
}

process sparcc_perm_cor{
  label 'sparcc'

  input:
  file perm from PERMS.flatten()

  output:
  file "perm_cor.txt" into PERMCORS

  """
  SparCC.py $perm -i ${params.iter} --cor_file=perm_cor.txt
  """
}

process sparcc_pval{
  label 'sparcc'
  publishDir "${params.outdir}/pvals", mode: 'rellink'

  input:
  file "perm_cor_*.txt" from PERMCORS.collect()
  file cor from COR

  output:
  file 'pvals.txt'

  """
  PseudoPvals.py $cor perm_cor_#.txt ${params.perms} -o pvals.txt
  """
}

/*
process{
  executor = 'slurm'
  withLabel: 'sparcc' {
    module = 'anaconda'
    conda = "/opt/modules/pkgs/anaconda/3.6/envs/sparcc"
    maxForks = 20
    env.PATH = "/home/sur/software/sparcc/yonatanf-sparcc-3aff6141c3f1/:$PATH"
  }
}
*/
