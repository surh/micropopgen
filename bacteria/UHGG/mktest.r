# (C) Copyright 2020 Sur Herrera Paredes
# 
# This file is part of micropopgen.
# 
# micropopgen is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# micropopgen is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with micropopgen.  If not, see <http://www.gnu.org/licenses/>.

# setwd("/cashew/users/sur/exp/fraserv/2020/today2/work/e4/e7352f79f8de68f5cf53cd1935a047/")
# setwd("/cashew/users/sur/exp/fraserv/2020/today2/work/f2/505baea646b2c3217bc0ae6730f25f/")
library(tidyverse)

#' Title
#'
#' @param eff_file 
#' @param feat_file 
#'
#' @return
#' @export
#'
#' @examples
read_snv_data <- function(eff_file, feat_file){
  snv_effs <- read_tsv(eff_file,
                       col_types = cols(ref_id = col_character(),
                                        ref_pos = col_number(),
                                        coding = col_number(),
                                        snp_effect = col_character(),
                                        substitution = col_character()))
  stop_for_problems(snv_effs)
  
  snv_feat <- read_tsv(feat_file,
                       col_names = c("ref_id", "ref_pos", "feat_type", "feat_id"),
                       col_types = cols(ref_id = col_character(),
                                        ref_pos = col_number(),
                                        feat_type = col_character(),
                                        feat_id = col_character()))
  stop_for_problems(snv_feat)
  
  if(nrow(snv_feat) > 0){
    # Need to do garbage collect to deal with some weird cons limit issue in R but not Rstudio
    gc()
    snv_effs %>%
      left_join(snv_feat, by = c("ref_id", "ref_pos")) %>%
      mutate(feat_type = replace(feat_type, is.na(feat_type), "IGR"))
  }else{
    snv_effs %>%
      mutate(feat_id = NA,
             feat_type = "IGR")
  }
}

#' Title
#'
#' @param vcf_file 
#'
#' @return
#' @export
#'
#' @examples
read_vcf_data <- function(vcf_file){
  read_tsv(vcf_file, comment = "##",
                  col_types = cols(`#CHROM` = col_character(),
                                   ID = col_number(),
                                   REF = col_character(),
                                   ALT = col_character(),
                                   QUAL = col_character(),
                                   FILTER = col_character(),
                                   INFO = col_character(),
                                   FORMAT = col_character(),
                                   .default = col_number()),
                  na = c("", "na", ".")) %>%
    select(-ID, -QUAL, -FILTER, -FORMAT) %>%
    rename(ref_id = "#CHROM",
           ref_pos = POS) %>%
    separate(INFO, sep = ";", into = c("N", "AF")) %>%
    mutate(N = str_remove(N, "NG=") %>% as.numeric,
           AF = str_remove(AF, "AF=") %>% as.numeric)
}

#' Title
#'
#' @param allele 
#' @param group 
#' @param min_size 
#'
#' @return
#' @export
#'
#' @examples
mkdist <- function(allele, group, min_size = 5){
  counts <- table(group)
  if(length(counts) != 2 || any(counts < 5)){
    return(NA)
  }
  
  tab <- table(allele, group)
  if(diag(tab) == 0 || sum(diag(tab)) == length(allele)){
    return("D")
  }else{
    return("P")
  }
}

#' Title
#'
#' @param dat 
#' @param genomes 
#' @param meta 
#' @param min_size 
#'
#' @return
#' @export
#'
#' @examples
test_mk <- function(dat, genomes, meta, min_size = 5){
  dat %>%
    # head(1000) %>%
    filter(feat_type == "CDS") %>%
    split(.$feat_id) %>%
    map_dfr(function(d, meta, genomes, min_size = 5){
      res <- d %>%
        select(ref_id, ref_pos, snp_effect, all_of(genomes)) %>%
        pivot_longer(cols = c(-ref_id, -ref_pos, -snp_effect),
                     names_to = "Genome", values_to = "allele") %>%
        filter(!is.na(allele)) %>%
        left_join(meta,
                  by = "Genome") %>%
        group_by(ref_pos) %>%
        summarise(snp_dist = mkdist(allele, group, min_size = min_size),
                  # ref_id = ref_id[1],
                  snp_effect = snp_effect[1],
                  .groups = 'drop')
      
      res <- table(factor(res$snp_dist, levels = c("D", "P")),
                   factor(res$snp_effect, levels = c("synonymous", "non-synonymous")))
      
      # res.test <- fisher.test(res)
      res <- as.vector(res)
      res.test <- fisher.test(matrix(c(res[3], res[1], res[4], res[2]), ncol = 2))
      
      tibble(Ps = res[2],
             Pn = res[4],
             Ds = res[1],
             Dn = res[3],
             OR = res.test$estimate,
             p.value = res.test$p.value)
    }, meta = meta,
    genomes = genomes,
    min_size = min_size, .id = "feat_id")
}

#' Title
#'
#' @param x 
#' @param m 
#'
#' @return
#' @export
#'
#' @examples
dummy_combn <- function(x, m = 2){
  if(length(x) >= m){
    combn(x, m = m)
  }else{
    matrix(ncol = 0, nrow = 1)
  }
}
#################################

opts <- commandArgs(trailingOnly = TRUE)
args <- list(snv_effects = opts[1],
             vcf = opts[2],
             snv_feats = opts[3],
             meta_file = opts[4],
             output = opts[5],
             min_size = as.numeric(opts[6]))

# args <- list(snv_effects = "MGYG-HGUT-00075_snveffs.tsv",
#              vcf = "MGYG-HGUT-00075.vcf.gz",
#              snv_feats = "MGYG-HGUT-00075_snvfeats.tsv",
#              meta_file = "metadata.txt",
#              output = "MGYG-HGUT-00075_mktest.tsv",
#              min_size = 5)
args

cat("Reading snv data...\n")
dat <-  read_snv_data(eff_file = args$snv_effects, feat_file = args$snv_feats)
cat("Reading genome ids...\n")
genomes <- (read_tsv(args$vcf, comment = "##",
                     n_max = 0,
                     col_types = cols(.default = col_character())) %>%
              colnames)[-(1:9)]
cat("Reading genotypes...\n")
dat <- dat %>%
  full_join(read_vcf_data(vcf_file = args$vcf),
            by = c("ref_id", "ref_pos"))

# ?HMVAR::midas_mktest
cat("Reading metadata..\n")
meta <- read_tsv(args$meta_file,
                 col_types = cols(CMseq = col_character())) %>%
  filter(Genome %in% genomes) %>%
  select(Genome, Genome_type, Country, Continent) %>%
  filter(!is.na(Continent))
genomes <- meta$Genome
# meta
continents <- unique(meta$Continent)
# continents
table(meta$Continent, useNA = "always")

# mktest <- tibble()
if(length(continents) == 1){
  cat("Not enough continents\n")
}else if(length(continents) == 2){
  cat("Two continents..\n")
  if(all(table(meta$Continent) >= args$min_size)){
    mktest <- test_mk(dat = dat, genomes = genomes,
                      meta = meta %>% select(Genome, group = Continent),
                      min_size = args$min_size) %>%
      mutate(group1 = continents[1], group2 = continents[2])
    # mktest
    mktest %>%
      filter(Ps + Pn + Dn + Ds > 0) %>%
      write_tsv(args$output)
  }else{
    cat("\tContinents didn't pass threshold...\n")
  }
}else if(length(continents) > 2){
  cat("More than two continents...\n")
  mktest <- bind_rows(continents %>%
                        map_dfr(function(continent, dat, meta, genomes, min_size = 5){
                          if(sum(meta$Continent == continent) < min_size){
                            return(NULL)
                          }
                          if(sum(meta$Continent != continent) < min_size){
                            return(NULL)
                          }
                          cat("Testing", continent, "vs other\n")
                          test_mk(dat = dat, genomes = genomes,
                                            meta = meta %>%
                                              select(Genome, group = Continent) %>%
                                              mutate(group = replace(group, group != continent, "other")),
                                            min_size = args$min_size) %>%
                            mutate(group1 = continent, group2 = "other")
                        }, dat = dat,
                        meta = meta,
                        genomes = genomes,
                        min_size = args$min_size),
                      table(meta$Continent) %>%
                        # For some reason as_tibble cannot convert directly the output of table
                        as.data.frame() %>%
                        as_tibble() %>%
                        filter(Freq >= args$min_size) %>%
                        select(Var1) %>%
                        unlist  %>%
                        as.character() %>%
                        dummy_combn(m = 2) %>%
                        t %>%
                        as_tibble(.name_repair = "minimal") %>%
                        pmap_dfr(function(V1, V2, dat, meta, genomes, min_size = 5){
                          cat("Testing", V1, "vs", V2, "\n")
                          test_mk(dat = dat,
                                  genomes = meta %>%
                                    filter(Continent %in% c(V1, V2)) %>%
                                    select(Genome) %>% unlist %>% as.character(),
                                  meta = meta %>%
                                    select(Genome, group = Continent) %>%
                                    filter(group %in% c(V1, V2)),
                                  min_size = min_size) %>%
                            mutate(group1 = V1, group2 = V2)
                        }, dat = dat,
                        meta = meta,
                        genomes = genomes,
                        min_size = args$min_size))
  
  # Only write if there is output
  if(nrow(mktest) > 0){
    mktest %>%
      filter(Ps + Pn + Dn + Ds > 0) %>%
      write_tsv(args$output)
  }else{
    cat("\tNo valid comparisons\n")
  }
}

# # mktest
# 
# if(nrow(mktest) > 0){
# 
# }
# 
# 
# mktest %>%
#   filter(Ps + Pn + Dn + Ds > 0) %>%
#   arrange(p.value) %>%
#   ggplot(aes(x = p.value)) +
#   geom_histogram(bins = 20)
# 
# mktest %>%
#   filter(Ps + Pn + Dn + Ds > 0) %>%
#   arrange(p.value) %>%
#   print(n = 30)
