library(ggplot2)
library(reshape2)
library(GGally)

Tab <- read.table("Porphyromonas_sp_57899_mk_results.txt", sep = "\t", header = TRUE)
head(Tab)

# Number of genes
nrow(Tab)

# Plot histograms of all p-values
Pvals <- Tab[,c("gene","hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
Pvals <- melt(Pvals, id.vars = "gene", value.name = "p.value", variable.name = "Test")
p1 <- ggplot(Pvals, aes(x = p.value)) +
  facet_wrap(~ Test, nrow = 2) +
  geom_histogram(bins = 20) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
ggsave("all_pvals.hist.svg", p1, width = 8, height = 4)

# Plot only those with at least one fixed non-synonymous mutation
Pvals <- subset(Tab, Dn > 0)
nrow(Pvals)
Pvals <- Pvals[,c("gene","hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
Pvals <- melt(Pvals, id.vars = "gene", value.name = "p.value", variable.name = "Test")
p1 <- ggplot(Pvals, aes(x = p.value)) +
  facet_wrap(~ Test, nrow = 2) +
  geom_histogram(bins = 20) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
ggsave("withDn_pvals.hist.svg", p1, width = 8, height = 4)


### Barplot number of tests
Pvals <- Tab[,c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
Pvals <- apply(is.na(Pvals),2,function(x) table(factor(x,levels=c(TRUE,FALSE))))
Pvals <- melt(Pvals, value.name = "Count", varnames = c("is.na","test"))
p1 <- ggplot(Pvals, aes(x = test, y = Count)) + 
  geom_bar(aes(fill = is.na), stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90))
p1
ggsave("sucess_tests.barplot.svg", p1, width = 6, height = 4)

### Plot correlation all
#Pvals <- Tab[,c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
pval_columns <- c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")
pval_columns <- which(colnames(Tab) %in% pval_columns )
p1 <- ggpairs(subset(Tab, !is.na(g_none_p)),columns = pval_columns) + theme(axis.text.x = element_text(angle = 90))
p1
ggsave("correlation_pvals_all.svg", p1, width = 10, height = 10)

### Plot correlation with Dn
Pvals <- subset(Tab, Dn > 0)
#Pvals <- Pvals[,c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
pval_columns <- c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")
pval_columns <- which(colnames(Pvals) %in% pval_columns )
p1 <- ggpairs(subset(Pvals, !is.na(g_none_p)),columns = pval_columns) + theme(axis.text.x = element_text(angle = 90))
p1
ggsave("correlation_pvals_withDn.svg", p1, width = 10, height = 10)

## Plot two versions of neutrality index
# head(Tab)
# ggplot(Tab, aes(x = ni, y = ratio))


p1 <- ggplot(Tab, aes(x = ni, y = -log10(g_none_p))) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_point()
p1


