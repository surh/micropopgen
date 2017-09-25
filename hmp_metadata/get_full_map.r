tab <- read.table("~/micropopgen/data/HMP/hmp_catalogue_wgs_all_metadata.csv", header = TRUE, sep = ",")
head(tab)
colnames(tab)

tab <- tab[,c("Sequence.Read.Archive.ID", "HMP.Isolation.Body.Site","HMP.Isolation.Body.Subsite")]
head(tab)
tab$HMP.Isolation.Body.Subsite[ tab$HMP.Isolation.Body.Site == "gastrointestinal_tract" ] <- "Stool"

ftable(HMP.Isolation.Body.Site ~ HMP.Isolation.Body.Subsite, tab)
table(tab$HMP.Isolation.Body.Subsite)

tab <- tab[,c(1,3)]
tab <- tab[ tab$HMP.Isolation.Body.Subsite != "", ]
colnames(tab) <- c("ID","Group")
head(tab)

subset(tab, ID %in% names(which(table(tab$ID) > 1)))

tab <- lapply(unique(tab$ID), function(x){
  temp <- subset(tab, ID %in% x)
  if(nrow(temp) > 1)
    temp <- temp[1,]
  
  temp
})
tab <- do.call(rbind,tab)

head(tab)

write.table(tab,"full_map.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
