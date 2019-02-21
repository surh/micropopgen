library(tidyverse)

args <- list(file = "checkm_results.txt")


Res <- read_tsv(file = args$file, col_types = 'ccnnnnnnnnnnnn')


p1 <- ggplot(Res, aes_string(x = "Completeness")) +
  geom_density() +
  scale_x_log10()
p1

breaks <- c(0, 80, 90,95,96,97,98,99,100)
p1 <- ggplot(tibble(Completeness = cut(Res$Completeness,
                                       breaks,
                                       include.lowest = TRUE)),
             aes(x = Completeness)) +
  geom_histogram(stat = "count") +
  theme_classic()
p1

p1 <- ggplot(tibble(Contamination = cut(Res$Contamination,
                                       100 - breaks,
                                       include.lowest = TRUE, right = FALSE)),
             aes(x = Contamination)) +
  geom_histogram(stat = "count") +
  theme_classic()
p1

Res %>% filter(Contamination > 100)


breaks <- c(0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
p1 <- ggplot(tibble(Strain.heterogeneity = cut(Res$`Strain heterogeneity`,
                                               breaks,
                                               include.lowest = TRUE)),
             aes(x = Strain.heterogeneity)) +
  geom_histogram(stat = "count") +
  theme_classic()
p1




