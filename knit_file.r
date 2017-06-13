#!/usr/bin/env Rscript
library(knitr)

rmdfile <- commandArgs(trailingOnly = TRUE)[1]
opts <- commandArgs(trailingOnly = TRUE)[-1]

# Set working directory to directory where file was called
opts_knit$set(root.dir = getwd())
#knit2html(input = rmdfile)

rmarkdown::render(input = rmdfile,output_dir = getwd())
