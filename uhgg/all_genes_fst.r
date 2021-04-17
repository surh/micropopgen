#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
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

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Calculate Fst for all sites and features between",
                        "countries for UHGG SNV catalogue."))

  # Positional arguments
  p <- add_argument(p, "snvs",
                    help = paste("SNVs file from UHGG catalogue"),
                    type = "character")
  p <- add_argument(p, "map",
                    help = "Metadata file")

  # Optional arguments
  p <- add_argument(p, "--gff",
                     help = paste("Gene feature file (GFF3) for the genome"),
                     type = "character",
                     default = NULL)
  p <- add_argument(p, "--pop_col",
                    help = "Column with populations for Fst",
                    type = "character",
                    default = "Country")
  p <- add_argument(p, "--pops",
                    help = "Populations to include in Fst calculations",
                    nargs = 2,
                    default = c("United States", "China"))
  p <- add_argument(p, "--only_isolates",
                    help = "If passed only isolate genomes included.",
                    flag = TRUE)
  p <- add_argument(p, "--outdir",
                    help = "output directory.",
                    default = "output")
  p <- add_argument(p, "--min_genomes",
                    help = paste("Minimum number of genomes per population"),
                    type = "numeric",
                    default = 5)

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)

  return(args)
}

args <- process_arguments()
# args <- list(snvs = "~/micropopgen/exp/2021/today/snv_catalogue/MGYG-HGUT-00022_snvs.tsv",
#              map = "~/micropopgen/exp/2021/today/test_meta.tsv",
#              gff = "~/micropopgen/exp/2021/today/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00022/genome/MGYG-HGUT-00022.gff",
#              pop_col = "Country",
#              pops = c("United States", "China"),
#              only_isolates = FALSE,
#              outdir = "output/",
#              min_genomes = 5)
print(args)
# q()
library(tidyverse)
library(HMVAR)

# Read genome metadata and select genomes from specified populations
cat("Reading genome metadata...\n")
Meta <- read_tsv(args$map,
                 col_types = cols(CMseq = col_number()))
Meta <- Meta %>%
  filter(.data[[args$pop_col]] %in% args$pops ) %>%
  filter(!is.na(.data[[args$pop_col]]))
if(args$only_isolates){
  cat("Cat selecting isolates...\n")
  Meta <- Meta %>%
    filter(Genome_type == "Isolate")
}

# Read SNVs
cat("Reading SNVs...\n")
snvs <- read_tsv(args$snvs)
genomes <- colnames(snvs)
genomes <- setdiff(genomes, c("Contig", "Pos", "Ref", "Alt"))

# Select genomes in SNV table
cat("Selecting final set of genomes...\n")
Meta <- Meta %>%
  filter(Genome %in% genomes)
genomes <- Meta$Genome

if(length(genomes) == 0){
  cat("There are no isolate genomes...\n")
  q()
}

# Count genomes per group
cat("Counting genomes per group...\n")
n_genomes <- Meta %>%
  group_by(.data[[args$pop_col]]) %>%
  summarise(n = length(MGnify_accession))

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Write number of genomes
cat("Writing genomes per group...\n")
filename <- file.path(args$outdir, "n_genomes.tsv")
write_tsv(n_genomes, filename)

# Checking if all groups have enough genomes
if(!all(n_genomes$n >= args$min_genomes)){
  cat("\tNot enough genomes in all groups. Finishing\n")
  q()
}
if(length(n$genomes_n) < 2){
  cat("\tNot enough groups for Fst calculation. Finishing\n")
  q()
}

# Select genomes from the right pops in snvs table and convert missing values
snvs <- snvs %>%
  select(Contig, Pos, Ref, Alt, all_of(genomes)) %>%
  mutate_at(.vars = vars(-Contig, -Pos, -Ref, -Alt),
            .funs = function(x){ifelse(x == 255, NA, x)})

cat("Preparing data for HMVAR...\n")
# Creating midas dat object
# freq <- snvs[1:1000,]
freq <- snvs

# Is this needed?
# if(nrow(freq) == 0){
#   cat("No SNVs in genomes")
#   q()
# }

info <- freq %>%
  transmute(site_id = paste0(Contig, ".", Pos),
            ref_id = Contig,
            ref_pos = Pos)
freq <- freq %>%
  mutate(site_id = paste0(Contig, ".", Pos)) %>%
  select(site_id, everything()) %>%
  select(-Contig, -Pos, -Ref, -Alt)
depth <- freq %>%
  mutate_at(vars(-site_id),
            .funs = function(x){
              replace(x, !is.na(x), 1)}) %>%
  mutate_at(.vars = vars(-site_id), .funs = function(x){
    replace_na(x, 0)})

# Site-level Fst
cat("Calculating site-level Fst...\n")
fst <- calculate_fst(Dat = list(freq = freq, depth = depth, info = info),
                     map = Meta %>%
                       select(sample = Genome, Group = .data[[args$pop_col]]),
                     support_thres = 1, verbose = FALSE)
cat("Writing site-level Fst...\n")
filename <- file.path(args$outdir, "site_fst.tsv")
write_tsv(fst$fst, filename)

if(!is.na(args$gff)){
  # Read GFF
  gff <- read_tsv(args$gff,
                  col_names = c("contig", "source", "feature",
                                "start", "end", "score",
                                "strand", "phase", "attributes")) %>%
    mutate(ID = attributes %>%
             map_chr(function(att){
               ((att %>%
                   str_split(";"))[[1]])[1] %>%
                 str_remove("^ID=")
             }))

  cat("Calculating feature level fst...\n")
  gff <- gff %>%
    select(gene_id = ID, contig, start, end, strand, feature) %>%
    mutate(Fst = NA,
           n_sites = NA)
  for(f in 1:nrow(gff)){
    dat <- fst$fst %>%
      filter(ref_id == gff$contig[f]) %>%
      filter(ref_pos >= gff$start[f] & ref_pos <= gff$end[f]) %>%
      filter(!is.na(Fst))

    if(nrow(dat) > 0){
      gff$Fst[f] <- max(0, sum(dat$a)/sum(dat$a + dat$b + dat$c))
      gff$n_sites[f] <- nrow(dat)
    }else{
      gff$n_sites[f] <- 0
    }
  }
  cat("Writing feature Fst")
  filename <- file.path(args$outdir, "feature_fst.tsv")
  write_tsv(gff, filename)
}

#
#   # Read SNVs and select genomes
#   # snvs <- read_tsv(snv_file)
#   # snvs <- snvs %>%
#   #   select(Contig, Pos, Ref, Alt, all_of(intersect(meta$Genome, colnames(snvs))))
#
#
#   # # Read RERperms and select genes
#   # rerperms <- read_tsv(f)
#   # set.seed(6523)
#   # rer_genes <- rerperms %>%
#   #   filter(!is.na(FDR)) %>%
#   #   filter(FDR < args$fdr_thres) %>%
#   #   select(gene_id) %>%
#   #   unlist %>%
#   #   as.character()
#   # baseline_genes <- rerperms %>%
#   #   filter(!is.na(FDR)) %>%
#   #   filter(FDR > 0.2) %>%
#   #   select(gene_id) %>%
#   #   unlist %>%
#   #   as.character()
#   # baseline_genes <- sample(baseline_genes, size = length(rer_genes), replace = FALSE)
#   #
#   # n_genomes <- n_genomes %>%
#   #   bind_rows(tibble(spec = spec,
#   #                    n_USA = sum(meta$Country == "USA"),
#   #                    n_CHN = sum(meta$Country == "CHN"),
#   #                    rer_genes = length(rer_genes),
#   #                    baseline_genes = length(baseline_genes)))
#
#
#
#   # Get only CDS
#   gff <- gff %>%
#     filter(feature == "CDS") %>%
#     select(contig, start, end, ID)
#
#   # Fst
#   Res <- gff %>%
#     rename(gene_id = ID) %>%
#     filter(gene_id %in% c(rer_genes, baseline_genes)) %>%
#     mutate(g_Fst = NA)
#   # Res
#   cat("\tCalculating Fst\n")
#   for(g in 1:nrow(Res)){
#     # g <- 1
#
#
#   }
#   # Join output and write
#   Res <- Res %>%
#     left_join(rerperms %>%
#                 select(gene_id, Rho, fdr_rer = FDR), by = "gene_id")
#
#   filename <- file.path("uhgg_fst_rer/", paste0(spec, ".tsv"))
#   write_tsv(Res, filename)
# }
#
# write_tsv(n_genomes, "uhgg_n_genomes.tsv")
# # p1 <- ggplot(Res, aes(x = Rho, y = g_Fst)) +
# #   geom_point(aes(col = fdr_rer < 0.05)) +
# #   theme_classic()
# # p1
# #
# # p1 <- ggplot(Res, aes(x = fdr_rer < 0.05, y = g_Fst)) +
# #   geom_violin(draw_quantiles = c(0.75, 0.5, 0.25)) +
# #   theme_classic()
# # p1
