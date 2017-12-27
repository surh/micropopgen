library(ggplot2)

set.seed(133)
dat <- rbind(data.frame(Bacteria = "No Bacteria",
                        Genotype = c("WT_h","WT_p","WT_l","h->p","p->h","l->p","l->h"),
                        Fitness = c(5,5.1,4.8,5.2,4.9,5,5.1),
                        Background = c("Human","Plant","Lab","Human","Plant","Lab","Lab")),
             data.frame(Bacteria = "SynCom_h",
                        Genotype = c("WT_h","WT_p","WT_l","h->p","p->h","l->p","l->h"),
                        Fitness = c(10,3,2,6,5.5,2.4,5),
                        Background = c("Human","Plant","Lab","Human","Plant","Lab","Lab")),
             data.frame(Bacteria = "SynCom_p",
                        Genotype = c("WT_h","WT_p","WT_l","h->p","p->h","l->p","l->h"),
                        Fitness = c(3.1,9.8,3,5.7,6,5.5,2.8),
                        Background = c("Human","Plant","Lab","Human","Plant","Lab","Lab"))
             
)
# dat$Genotype <- factor(dat$Genotype, levels = c("WT_h","WT_p","WT_l","h->p","p->h","l->p","l->h"))
dat$Genotype <- factor(dat$Genotype, levels = c("WT_h","h->p","WT_p","p->h","WT_l","l->p","l->h"))

dat

p1 <- ggplot(dat,aes(y = Fitness, x = Genotype)) +
  facet_grid(~ Bacteria) +
  geom_bar(aes(fill = Background), stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, color = "black"))
p1
ggsave("~/micropopgen/docs/proposals/2017_HHMI_HannaHGray/crispey_cartoon.svg",p1,width = 8, height = 4)
