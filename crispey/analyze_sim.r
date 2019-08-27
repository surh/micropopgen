library(tidyverse)
library(ggplot2)
library(lme4)
library(AMOR)


# Simulate and save
evo <- tournament_selection(pop = pop, k = k, G = max(timepoints))
# save(evo, file = "sim_pops.rdat")
counts <- pops_count_table(evo)
# save(counts, file = "sim_counts.rdat")

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

