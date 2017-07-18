v13 <- read.table("~/micropopgen/data/HMP/ppAll_V13_map.txt", sep = "\t", header = TRUE)
head(v13)

v35 <- read.table("~/micropopgen/data/HMP/ppAll_V35_map.txt", sep = "\t", header = TRUE)
head(v35)

v69 <- read.table("~/micropopgen/data/HMP/ppAll_V69_map.txt", sep = "\t", header = TRUE)
head(v69)

dat <- rbind(v13,v35,v69)

Res <- NULL
for(sample in levels(dat$SRS_SampleID)){
  # sample <- levels(dat$SRS_SampleID)[1]
  
  if (sample == "Unknown") next
  
  temp <- subset(dat,SRS_SampleID == sample)
  
  
  
  rsid <- unique(temp$RSID)
  sex <- unique(temp$Sex)
  body_subsite <- unique(temp$HMPBodySubsite)
  body_site <- unique(temp$HMPBodySite)
  visit_no <- unique(temp$VisitNo)
  
  res <- data.frame(SRS.ID = sample, RSID = rsid, Sex = sex, Body.Site = body_site, Body.Subsite = body_subsite, Visit.No = visit_no, stringsAsFactors = FALSE)
  if(nrow(res) > 1){
    stop("ERROR")
  }
  Res <- rbind(Res,res)
}
  
head(Res)  
write.table(Res, "subject_metadata.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

