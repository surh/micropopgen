library(tidyverse)
library(AMOR)

single_tournament <- function(parents, k = 2){
  ii <- sample(nrow(parents), size = k, replace = FALSE)
  tournament <- parents[ii, ]
  ii <- which.max(tournament$s)
  return(tournament[ii,])
}

tournament_round <- function(parents, k = 2, N = nrow(pop)){
  1:N %>%
    map_dfr(~single_tournament(parents = parents, k = k))
}

tournament_selection <- function(pop, k = 2, G = 1){
  N <- nrow(pop)
  
  Evo <- tibble(gen0 = pop$id)
  parents <- pop
  for(i in 1:G){
    cat("Generation: ", i, "\n")
    parents <- tournament_round(parents = parents, k = k, N = N)
    Evo <- Evo %>%
      bind_cols(!!paste0("gen", i) := parents$id)
  }
  
  return(Evo)
}

pops_count_table <- function(pops){
  # pops <- evo
  counts <- pops %>%
    map(function(x){
      tab <- table(x)
      tibble(id = names(tab), count = tab)
    }) %>%
    reduce(left_join, by = "id")
  names(counts) <- c("id", names(pops))
  counts <- counts %>% 
    replace_na(as.list(setNames(rep(0, ncol(counts) - 1), names(pops))))
  ids <- counts$id
  counts <- counts %>%
    dplyr::select(-id) %>%
    as.matrix()
  row.names(counts) <- ids
  
  return(counts)
}

##################################
set.seed(12345)
n_oligos <- as.numeric(opts[1])
mean_barcodes <- as.numeric(opts[2])
n_significant <- as.numeric(opts[3])
n_cells <- as.numeric(opts[4])
timepoints <- 20
k <- as.numeric(opts[5])
seed <- as.numeric(opts[6])

set.seed(seed)
cat(seed,"\n")

# Number of barcodes per oligo, add overdispersion?
n_barcodes <- rpois(n_oligos, lambda = mean_barcodes)
barcode_ids <- tibble(oligo = paste0('o.', rep(1:n_oligos, n_barcodes)),
                      barcode = paste0("b.",  1:sum(n_barcodes)))
barcode_ids$s <- 1
selected_oligos <- sample(unique(barcode_ids$oligo), size = n_significant, replace = FALSE)
selected_barcodes <- barcode_ids %>%
  filter(oligo %in% selected_oligos) %>%
  transmute(id = paste0(oligo, "_", barcode)) %>%
  unlist
barcode_ids$s[ barcode_ids$oligo %in% selected_oligos ] <- barcode_ids$s[ barcode_ids$oligo %in% selected_oligos ] + 1

# Create starting population
pop <- barcode_ids[ sample(x = 1:nrow(barcode_ids), size = n_cells, replace = TRUE), ] %>%
  mutate(id = paste0(oligo, "_", barcode)) %>%
  dplyr::select(id, s)

# Simulate and save
evo <- tournament_selection(pop = pop, k = k, G = max(timepoints))
save(evo, file = "sim_pops.rdat")
counts <- pops_count_table(evo)
save(counts, file = "sim_counts.rdat")
