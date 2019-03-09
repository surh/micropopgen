setwd("/godot/users/sur/exp/fraserv/2019/today/")
library(tidyverse)

args <- list(genome_dirs = "genome_dirs.txt",
             genome_meta = "patric_metadata/genome_metadata_notrailing",
             genome_lineage = "patric_metadata/genome_lineage")

# Read data
dirs <- read_tsv(args$genome_dirs, col_types = 'c', col_names = "path")
dirs

meta <- read_tsv(args$genome_meta,
                 col_types = cols(genome_id = col_character(),
                                  genome_name = col_character(),
                                  organism_name = col_character(),
                                  taxon_id = col_character(),
                                  strain = col_character(),
                                  bioproject_accession = col_character(),
                                  biosample_accession = col_character(),
                                  assembly_accession = col_character(),
                                  genbank_accessions = col_character(),
                                  refseq_accessions = col_character(),
                                  chromosomes = col_number(),
                                  plasmids = col_number(),
                                  contigs = col_number(),
                                  sequences = col_number(),
                                  genome_length = col_number(),
                                  gc_content = col_number(),
                                  patric_cds = col_number(),
                                  brc1_cds = col_number(),
                                  refseq_cds = col_number(),
                                  latitude = col_character(),
                                  longitude = col_character(),
                                  altitude = col_character(),
                                  depth = col_character(),
                                  host_age = col_character(),
                                  temperature_range = col_character(),
                                  optimal_temperature = col_character(),
                                  salinity = col_character(),
                                  .default = col_character()))
meta

lineage <- read_tsv(args$genome_lineage, col_types = cols(.default = col_character()))
lineage

# Process dirs
dirs <- dirs %>%
  mutate(genome_id = basename(.$path)) %>%
  mutate(dir = dirname(.$path)) %>%
  mutate(dir = basename(dir)) %>%
  mutate(dir = str_split(string = .$dir, pattern = "_", n = 2)) %>%
  rowwise() %>%
  mutate(dir_genus = unlist(dir)[1], dir_species = unlist(dir)[2]) %>%
  select(-dir)
dirs

table(dirs$dir_genus)

dirs <- dirs %>% left_join(lineage, by = "genome_id") %>%
  select(-taxon_lineage_ids, -taxon_lineage_names, -taxon_id)
dirs


dirs %>% filter(dir_genus != genus)
