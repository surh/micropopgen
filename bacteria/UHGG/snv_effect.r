# (C) Copyright 2020 Sur Herrera Paredes
# 
# This file is part of .
# 
#  is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
#  is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with .  If not, see <http://www.gnu.org/licenses/>.

# setwd("/cashew/users/sur/exp/fraserv/2020/today")
# setwd("/cashew/users/sur/exp/fraserv/2020/today/work/7c/01c0f6afbcf343b6ca6941768caddc")
# library(tidyverse)
library(magrittr)
library(PopGenome)

opts <- commandArgs(trailingOnly = TRUE)

vcf_dir <- opts[1]
gff_dir <- opts[2]
contig_fna <- opts[3]
output <- opts[4]

# vcf_dir <- "test/nosnvs/vcf/"
# gff_dir <- "test/nosnvs/gff/"
# contig_fna <- "test/nosnvs/fna/GUT_GENOME000001_100.fasta"

# vcf_dir <- "test/snvs/vcf/"
# gff_dir <- "test/snvs/gff/"
# contig_fna <- "test/snvs/fna/GUT_GENOME000001_1.fasta"

# vcf_dir <- "vcf/"
# gff_dir <- "gff/"
# contig_fna <- "GUT_GENOME000001_106.fasta"
# output <- "GUT_GENOME000001_106.tsv"

cat("========== params ==========\n")
cat(vcf_dir, "\n")
cat(gff_dir, "\n")
cat(contig_fna, "\n")
cat(output, "\n")
cat("========== params ==========\n")

vars <- readData(vcf_dir, format="VCF",
                 gffpath = gff_dir,
                 include.unknown = TRUE)

if(vars@n.biallelic.sites > 0){
  vars <- set.synnonsyn(vars, ref.chr = contig_fna)
  ctg <- vars@region.names %>% stringr::str_remove("[.]vcf$")
  
  snvs <- tibble::tibble(ref_id = ctg,
                         position = vars@region.data@biallelic.sites[[1]],
                         synonymous = vars@region.data@synonymous[[1]],
                         transitions = vars@region.data@transitions[[1]],
                         coding = vars@region.data@CodingSNPS[[1]]) %>%
    dplyr::mutate(snp_effect = replace(synonymous, synonymous == 1, "synonymous")) %>%
    dplyr::mutate(snp_effect = replace(snp_effect, synonymous == 0, "non-synonymous")) %>%
    dplyr::mutate(snp_effect = replace(snp_effect, is.na(synonymous), NA)) %>%
    dplyr::mutate(substitution = replace(transitions, transitions == 1, "transition")) %>%
    dplyr::mutate(substitution = replace(substitution, transitions == 0, "transversion")) %>%
    dplyr::mutate(substitution = replace(substitution,
                                         substitution != "transition" & substitution != "transversion",
                                         NA)) %>%
    dplyr::mutate(coding = 1*coding) %>%
    dplyr::select(ref_id, ref_pos = position, coding, snp_effect, substitution) 
  # snvs
  
  readr::write_tsv(snvs, path = output)
}

