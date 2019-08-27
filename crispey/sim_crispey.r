library(tidyverse)
library(ggplot2)
library(lme4)
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
n_oligos <- 100
mean_barcodes <- 10
n_significant <- 10
n_cells <- 1e4
timepoints <- c(0, 4, 8 , 12, 16, 20)
k <- 20

# n_oligos <- 100
# mean_barcodes <- 3
# n_significant <- 10
# n_cells <- 1e3
# timepoints <- c(0, 4, 8 , 12, 16, 20)
# k <- 20


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
save(evo, file = "sim_data.rdat")
counts <- pops_count_table(evo)
save(counts, file = "sim_counts.rdat")

# Create dataset
# Map <- data.frame(do.call(rbind, 
#                           row.names(counts) %>%
#                             map(~str_split(string = ., pattern = "_")[[1]])),
#                   stringsAsFactors = FALSE)
# row.names(map) <- row.names(counts)
Map <- data.frame(ID = colnames(counts),
                  time = as.numeric(str_replace(colnames(counts),
                                                pattern = "^gen", replacement = "")),
                  stringsAsFactors = FALSE)
row.names(Map) <- colnames(counts)
Dat <- create_dataset(Tab = counts, Map = Map)

# Plot
dat <- tibble()
for(b in intersect(selected_barcodes, taxa(Dat))){
  dat <- dat %>% bind_rows(plotgg_taxon(Dat, taxon = b, x = "time")$data %>%
                             bind_cols(barcode = rep(b, length(samples(Dat)))))
}
dat$oligo <- dat$barcode %>%
  str_split(pattern = "_") %>%
  map_chr(~.[1])
p1 <- ggplot(dat, aes(x = time, y = Abundance,
                      group = barcode)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_blackbox()
p1
p1 <- ggplot(dat, aes(x = time, y = Abundance,
                      group = barcode)) +
  facet_wrap(~ oligo) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_blackbox()
p1


######

evo %>%
  map(~as.matrix(table(.)))

set.seed(12345)
dat <- tibble(oligo = sample(1:n_oligos, size = N, replace = T),
              barcode = sample(1:n_barcodes, size = N, replace = T))
dat <- dat %>%
  transmute(oligo = paste0("o", oligo), barcode = paste0("b", barcode)) %>%
  mutate(id = paste0(oligo, ".", barcode))
dat

t0 <- table(dat$id) / N
Dat <- str_split(string = names(t0), pattern = "[.]", simplify = T) %>%
  as.tibble %>% 
  dplyr::select(oligo = V1, barcode = V2) %>%
  bind_cols(id = names(t0)) %>%
  mutate(s_coef = s_coef * 1 * (oligo %in% paste0("o", 1:n_sig))) %>%
  bind_cols(t0 = as.numeric(t0))
t_prev <- Dat$t0
for(t in 1:20){
  # t <- 1
  p <- t_prev + Dat$s_coef
  p <- N * p / sum(p)
  temp <- rpois(nrow(Dat), lambda = p)
  temp <- temp / sum(temp)
  if(t %in% timepoints){
    Dat[[paste0("t", t)]] <- temp
  }
  t_prev <- temp
}

Dat
# Dat <- Dat %>% filter(rowSums(Dat[,6:10] == 0) < 3)
Dat <- Dat %>% filter(t20 > 0)
Dat
Dat <- Dat %>% gather(timepoint, freq, t0:t20) %>%
  mutate(time = as.numeric(str_replace(timepoint,
                                       pattern = "^t",
                                       replacement = "")))
Dat
# p1 <- ggplot(subset(Dat, oligo == "o100"),
#              aes(x = time, y = freq, group = barcode)) +
#   geom_line(alpha = 0.2)
# p1

Dat$freq <- log2(Dat$freq * N + 1)
# m1 <- lmer(freq ~ time + oligo + (1|oligo) + (1|id) + (1|time), data = Dat, verbose = TRUE)
m1 <- lmer(freq ~ time + oligo + (1|oligo), data = Dat, verbose = TRUE)




p1 <- ggplot(subset(Dat, oligo == "o99"),
             aes(x = time, y = freq, group = barcode)) +
  geom_line(alpha = 0.2)
p1

p1 <- ggplot(subset(Dat, oligo == "o2"),
             aes(x = freq)) +
  facet_grid(~timepoint) +
  geom_density() +
  geom_vline(xintercept = 3)
p1

p1 <- ggplot(subset(Dat),
             aes(x = factor(timepoint, levels = c("t0","t4","t8","t12","t16","t20")),
                 y = freq, color = oligo %in% paste0("o", 1:n_sig))) +
  geom_boxplot()
p1


p1 <- ggplot(subset(Dat),
             aes(x = factor(timepoint, levels = c("t0","t4","t8","t12","t16","t20")),
                 y = freq, color = oligo %in% paste0("o", 1:n_sig))) +
  geom_violin()
p1
