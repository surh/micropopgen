setwd("/cashew/users/sur/exp/fraserv/2020/today")
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
  
  snv_feat <- read_tsv(feat_file,
                       col_names = c("ref_id", "ref_pos", "feat_type", "feat_id"),
                       col_types = cols(ref_id = col_character(),
                                        ref_pos = col_number(),
                                        feat_type = col_character(),
                                        feat_id = col_character()))
  
  
  snv_effs %>%
    left_join(snv_feat, by = c("ref_id", "ref_pos")) %>%
    mutate(feat_type = replace(feat_type, is.na(feat_type), "IGR"))
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
#################################

opts <- commandArgs(trailingOnly = TRUE)
args <- list(snv_effects = opts[1],
             vcf = opts[2],
             snv_feats = opts[3],
             meta_file = opts[4])

args <- list(snv_effects = "output/snvs/MGYG-HGUT-00001.tsv",
             vcf = "output/tabix/MGYG-HGUT-00001.vcf.gz",
             snv_feats = "output/snv_feats/MGYG-HGUT-00001.tsv",
             meta_file = "/cashew/shared_data/mgnify/v1.0/genomes-all_metadata.tsv")


dat <-  read_snv_data(eff_file = args$snv_effects, feat_file = args$snv_feats)
genomes <- (read_tsv(args$vcf, comment = "##",
                     n_max = 0,
                     col_types = cols(.default = col_character())) %>%
              colnames)[-(1:9)]
dat <- dat %>%
  full_join(read_vcf_data(vcf_file = args$vcf),
            by = c("ref_id", "ref_pos"))

# ?HMVAR::midas_mktest
meta <- read_tsv(args$meta_file,
                 col_types = cols(CMseq = col_character())) %>%
  filter(Genome %in% genomes) %>%
  select(Genome, Genome_type, Country, Continent)
meta


mkdist <- function(allele, group){
  tab <- table(allele, group)
  if(diag(tab) == 0 || sum(diag(tab)) == length(allele)){
    return("D")
  }else{
    return("P")
  }
}


mktest <- dat %>%
  filter(feat_type == "CDS") %>%
  split(.$feat_id) %>%
  map_dfr(function(d, meta, genomes){
    res <- d %>%
      select(ref_id, ref_pos, snp_effect, genomes) %>%
      pivot_longer(cols = c(-ref_id, -ref_pos, -snp_effect),
                   names_to = "Genome", values_to = "allele") %>%
      filter(!is.na(allele)) %>%
      left_join(meta %>%
                  select(Genome, group = Continent),
                by = "Genome") %>%
      group_by(ref_pos) %>%
      summarise(snp_dist = mkdist(allele, group),
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
  }, meta = meta, genomes = genomes, .id = "feat_id")
mktest
mktest %>%
  arrange(p.value) %>%
  print(n = 50)


