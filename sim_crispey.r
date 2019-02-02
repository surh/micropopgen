library(tidyverse)
library(ggplot2)
library(lme4)

n_oligos <- 100
n_barcodes <- 30 * n_oligos
n_sig <- 10
s_coef <- 1e-7
N <- 1000000
timepoints <- c(0, 4, 8 , 12, 16, 20)

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
