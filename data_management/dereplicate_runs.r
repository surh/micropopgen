files <- opts
# files <- c("SRA045646_map.txt", "SRA050230_map.txt")
print(files)

TAB <- NULL
for (f in files){
  t <- read.table(f, header = TRUE)
  TAB <- rbind(TAB,t)
}
head(TAB)
colnames(TAB) <- c("Run", "Sample")

nrows <- nrow(TAB)
nsamples <- length(table(TAB$Sample))
cat("There are ", nsamples, " samples, and ", nrows, " rows\n")

TAB <- as.data.frame(do.call(rbind, strsplit(x = unique(sort(paste(TAB$Run, TAB$Sample))), split = " ")))
colnames(TAB) <- c("Run", "Sample")
head(TAB)

nruns <- nrow(TAB)
cat("There are ", nsamples, " samples, and ", nruns, " runs\n")

write.table(TAB, "runs_map.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)