library(tidyverse)
library(ggplot2)
library(lme4)

single_tournament <- function(parents, k = 2){
  ii <- sample(nrow(parents), size = k, replace = FALSE)
  return((parents[ii, ] %>% arrange(desc(s)))[1,])
}

tournament_round <- function(parents, k = 2, N = nrow(pop)){
  1:N %>%
    map_dfr(~single_tournament(parents = parents, k = k))
  }

tournament_selection <- function(pop, k = 2, N = nrow(pop), G = 1){
  tibble(gen0 = pop$id) %>% 
    bind_cols(setNames(1:G, paste0("gen", 1:G)) %>%
                map(~tournament_round(parents = pop, k = k, N = N)) %>%
                map(~.$id) %>%
                bind_cols())
  
}

set.seed(12345)
n_oligos <- 1000
mean_barcodes <- 30
n_significant <- 10
n_cells <- 1e6
timepoints <- c(0, 4, 8 , 12, 16, 20)
k <- 100


# Number of barcodes per oligo, add overdispersion?
n_barcodes <- rpois(n_oligos, lambda = mean_barcodes)
# barcode_ids <- paste0(paste0("o.", rep(1:n_oligos, n_barcodes)),'_b.', 1:sum(n_barcodes))
barcode_ids <- tibble(oligo = paste0('o.', rep(1:n_oligos, n_barcodes)),
                      barcode = paste0("b.",  1:sum(n_barcodes)))
barcode_ids$s <- 1
selected_oligos <- sample(unique(barcode_ids$oligo), size = n_significant, replace = FALSE)
# selected_oligos
barcode_ids$s[ barcode_ids$oligo %in% selected_oligos ] <- barcode_ids$s[ barcode_ids$oligo %in% selected_oligos ] + 1
# barcode_ids %>%
#   print(n=100)

# Create starting population
pop <- barcode_ids[ sample(x = 1:nrow(barcode_ids), size = n_cells, replace = TRUE), ] %>%
  mutate(id = paste0(oligo, "_", barcode)) %>%
  select(id, s)

evo <- tournament_selection(pop = pop, k = k, N = n_cells, G = max(timepoints))
evo_counts <- evo %>%
  map(function(x){
    tab <- table(x)
    tibble(id = names(tab), count = tab)
  }) %>%
  reduce(left_join, by = "id")
names(evo_counts) <- c("id", names(evo))
evo_counts


selected_oligos
oligo_barcodes <- barcode_ids %>%
  filter(oligo == "o.8") %>%
  transmute(id = paste0(oligo, "_", barcode)) %>%
  unlist
evo_counts %>%
  filter(id %in% oligo_barcodes)
pop %>%
  filter(id %in% oligo_barcodes)


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
  select(oligo = V1, barcode = V2) %>%
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
