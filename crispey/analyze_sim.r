library(tidyverse)
library(ggplot2)
library(lme4)
library(AMOR)

counts_file <- "oligos1000_barcodes30_sig10_cells10000_k2_rep1/sim_counts.rdat"
barcodes_file <- "oligos1000_barcodes30_sig10_cells10000_k2_rep1/barcode_ids.rdat"
timepoints <- c(0, 4, 6, 12, 20)

load(file = counts_file)

# Create Dataset
Map <- data.frame(ID = colnames(counts),
                  time = as.numeric(str_replace(colnames(counts),
                                                pattern = "^gen", replacement = "")),
                  stringsAsFactors = FALSE)
row.names(Map) <- colnames(counts)
Dat <- create_dataset(Tab = counts, Map = Map)

# Find selected ids
load(barcodes_file)
selected_barcodes <- barcode_ids %>%
  dplyr::filter(s > 1) %>%
  pmap_chr(function(oligo, barcode, ...){ paste0(oligo, "_", barcode)})

# barcode_ids %>%
#   dplyr::filter(s > 1)
# barcode_ids %>%
#   dplyr::filter(s > 1) %>%
#   dplyr::select(oligo) %>%
#   table

colSums(remove_taxons(Dat, setdiff(taxa(Dat), selected_barcodes))$Tab)

# Plot
dat <- tibble()
for(b in intersect(selected_barcodes, taxa(Dat))){
  dat <- dat %>% bind_rows(plotgg_taxon(Dat, taxon = b, x = "time")$data %>%
                             bind_cols(barcode = rep(b, length(samples(Dat)))))
}
dat$oligo <- dat$barcode %>%
  str_split(pattern = "_") %>%
  map_chr(~.[1])
# p1 <- ggplot(dat, aes(x = time, y = Abundance,
#                       group = barcode)) +
#   geom_point() +
#   geom_smooth(method = "loess") +
#   theme_blackbox()
# p1
p1 <- ggplot(dat, aes(x = time, y = Abundance,
                      group = barcode)) +
  facet_wrap(~ oligo, scales = "free_y") +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_blackbox()
p1

# Select timepoints
Dat <- subset(Dat, time %in% timepoints, clean = TRUE)
# Select oligos with multiple timepoints
Dat <- remove_taxons(Dat, taxa(Dat)[ rowSums(Dat$Tab > 0) <= 3])
Dat

dat <- as_tibble(Dat$Tab) %>%
  bind_cols(barcode = taxa(Dat)) %>%
  gather(key = generation, value = count, -barcode) %>%
  left_join(Dat$Map %>%
              dplyr::select(generation = ID, time),
            by = "generation") %>%
  mutate(oligo = str_split(string = barcode, "_") %>% map_chr(~.[1]))
dat

m1 <- lmer(count ~ time + (1 | barcode) + (1 | oligo), data = dat %>% filter(barcode %in% selected_barcodes))
m1 <- lmer(count ~ time + (1 | barcode) + (1 | oligo), data = dat)
m1 <- lmer(count ~ time + (barcode | oligo), data = dat)
summary(m1)

lattice::dotplot(ranef(m1, which = "oligo"))
lattice::dotplot(ranef(m1))

o <- unique(dat$oligo)[2]
m1 <- lmerTest::lmer(count ~ time + (1|barcode) + (1|time), dat %>% filter(oligo == o))
summary(m1)
anova(m1)
lmerTest::ranova(m1)
