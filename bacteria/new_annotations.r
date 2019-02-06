library(tidyverse)

genes_file <- "hmp.sig.genes.GO.txt"
annot_dir <- "/godot/users/sur/exp/fraserv/2018/today4/CDS/"

genes <- read_tsv(genes_file)
genes

gene_annotations <- NULL
for(i in 1:nrow(genes)){
  # i <- 10
  
  spec <- genes$name[i]
  filename <- paste0(annot_dir, "/", spec,  ".emapper.annotations")
  
  if(file.exists(filename)){
    cat("\t",filename, "\n")
    
    annot <- read_tsv(filename, col_types = 'ccdnccccccccc', na = 'NA', comment = "#",
                      col_names = c("query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue", "seed_ortholog_score",
                                    "predicted_gene_name", "GO_terms", "KEGG_KOs", "BiGG_reactions", "Annotation_tax_scope",
                                    "OGs", "bestOG.evalue.score", "COG.cat", "eggNOG.annot"))
    annot <- annot %>% mutate(query_name = str_replace(query_name, "\\([\\+\\-]\\)_[0-9]+", ""))
    annot
    
    gene_id <- genes$gene[i]
    
    annot %>% filter(query_name == gene_id) %>% select()
    gene <- annot %>% filter(query_name == gene_id) %>%
      select(gene = query_name, predicted_gene_name, GO_terms, KEGG_KOs, BiGG_reactions, OGs, COG.cat, eggNOG.annot)
    
    gene_annotations <- gene_annotations %>% bind_rows(gene)
  }
}
genes <- genes %>% left_join(gene_annotations, by = "gene")
write_tsv(genes, "hmp_new_incomplete_annotations.txt")



rm(list=ls())

##################

genes_file <- "qin2012.sig.genes.GO.txt"
annot_dir <- "/godot/users/sur/exp/fraserv/2018/today4/CDS/"

genes <- read_tsv(genes_file)
genes

gene_annotations <- NULL
for(i in 1:nrow(genes)){
  # i <- 10
  
  spec <- genes$name[i]
  filename <- paste0(annot_dir, "/", spec,  ".emapper.annotations")
  
  if(file.exists(filename)){
    cat("\t",filename, "\n")
    
    annot <- read_tsv(filename, col_types = 'ccdnccccccccc', na = 'NA', comment = "#",
                      col_names = c("query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue", "seed_ortholog_score",
                                    "predicted_gene_name", "GO_terms", "KEGG_KOs", "BiGG_reactions", "Annotation_tax_scope",
                                    "OGs", "bestOG.evalue.score", "COG.cat", "eggNOG.annot"))
    annot <- annot %>% mutate(query_name = str_replace(query_name, "\\([\\+\\-]\\)_[0-9]+", ""))
    annot
    
    gene_id <- genes$gene[i]
    
    annot %>% filter(query_name == gene_id) %>% select()
    gene <- annot %>% filter(query_name == gene_id) %>%
      select(gene = query_name, predicted_gene_name, GO_terms, KEGG_KOs, BiGG_reactions, OGs, COG.cat, eggNOG.annot)
    
    gene_annotations <- gene_annotations %>% bind_rows(gene)
  }
}
genes <- genes %>% left_join(gene_annotations, by = "gene")
write_tsv(genes, "qin2012_new_incomplete_annotations.txt")
