// Params
params.midas_db = ''
params.bed_dir = ''
params.genomes_file = ''
params.outdir = 'output/'

// process params
genomes_file = file(params.genomes_file)

reader = genomes_file.newReader()
GENOMES = []
while( line = reader.readLine() ) {
  GENOMES = GENOMES + [tuple(line,
    file("${params.bed_dir}/${line}.bed"),
    file("${params.midas_db}/rep_genomes/${line}/genome.features.gz"))]
}

process bed_intersect_genes{
  module "anaconda"
  conda "/opt/modules/pkgs/anaconda/3.6/envs/fraserconda"

  publishDir params.outdir, mode: 'rellink', pattern: "*.selected.bed"
  stageInMode 'rellink'
  maxForks 100

  input:
  set genome, bed_file, feat_file from GENOMES

  output:
  file 'genome.features.bed'
  file "${genome}.selected.bed" optional true

  script:
  """
  # Convert features to BED
  zcat ${feat_file} | awk '{print \$2 "\t" \$3 "\t" \$4 "\t" \$1}' | \
    grep -v scaffold_id | sort -k1,1 -k2,2n > genome.features.bed

  # Find closest
  bedtools intersect -wa -a genome.features.bed \
    -b ${bed_file} > ${genome}.selected.bed
  """
}
