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
params.input = ''
params.iter = 20
params.perms = 100
params.outdir = 'output'

input = file(params.input)

process sparcc_cor{
  label 'sparcc'
  label 'highmem'
  publishDir "${params.outdir}/cor", mode: 'rellink'

  input:
  file input

  output:
  file 'cor_mat_SparCC.out' into COR
  file 'cov_mat_SparCC.out'

  """
  python ${workflow.projectDir}/SparCC.py $input -i ${params.iter}
  """
}

process sparcc_bootstraps{
  label 'sparcc'
  label 'standard'
  publishDir "${params.outdir}", mode: 'rellink'

  input:
  file input

  output:
  file 'perms/permutation_*.txt' into PERMS

  """
  python ${workflow.projectDir}/MakeBootstraps.py $input -n ${params.perms} -t permutation_#.txt -p perms/
  """
}

process sparcc_perm_cor{
  label 'sparcc'
  label 'highmem'

  input:
  file perm from PERMS.flatten()

  output:
  file "perm_cor.txt" into PERMCORS

  """
  python ${workflow.projectDir}/SparCC.py $perm -i ${params.iter} --cor_file=perm_cor.txt
  """
}

process sparcc_pval{
  label 'sparcc'
  label 'highmem'
  publishDir "${params.outdir}/pvals", mode: 'rellink'

  input:
  file "perm_cor_*.txt" from PERMCORS.collect()
  file cor from COR

  output:
  file 'pvals.txt'

  """
  # Need to start the permutations with 0-index.
  ln -s perm_cor_${params.perms}.txt perm_cor_0.txt
  python ${workflow.projectDir}/PseudoPvals.py $cor perm_cor_#.txt ${params.perms} -o pvals.txt
  """
}

// Example nexflow.config
/*
process{
  withLabel: 'sparcc' {
    module = 'anaconda'
    conda = "/opt/modules/pkgs/anaconda/3.6/envs/sparcc"
    maxForks = 300
  }
  withLabel: 'highmem'{
    memory = '20G'
    time = '2d'
  }
  withLabel: 'standard'{
    memory = '1G'
    time = '1d'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
}
*/
