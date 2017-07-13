


file <- "~/micropopgen/data/HMP/HMIWGS_healthy_run_metadata/all_runs.txt"
all_runs <- read.table(file,header = TRUE,sep = "\t")
head(all_runs)

file <- "~/micropopgen/data/HMP/HMIWGS_healthy.csv"
hmiwgs_samples <- read.table(file,sep = ",", quote = '"', header = TRUE)
head(hmiwgs_samples)

file <- "~/micropopgen/data/HMP/HMASM.csv"
hmasm_samples <- read.table(file,sep = ",", quote = '"', header = TRUE)
head(hmasm_samples)

file <- "~/micropopgen/data/HMP/HMASM-690.csv"
hmasm_samples_690 <- read.table(file,header = TRUE,sep = ",",quote = '"')
head(hmasm_samples_690)

# Compare body sisetes from hmasm sets
sites1 <- as.character(hmasm_samples_690$Body.Site)
names(sites1) <- as.character(hmasm_samples_690$SRS.ID)
sites2 <- as.character(hmasm_samples$Body.Site)
names(sites2) <- as.character(hmasm_samples$SRS.ID)

head(sites1)
head(sites2)

shared <- intersect(names(sites1),names(sites2))

if(!all(sites1[ shared ] == sites2[ shared ])){
  stop("ERROR: Not all samples match on body site")
}else{
  cat("All samples match\n")
}

hmasm_samples$HMASM.QC <- !is.na(sites1[ as.character(hmasm_samples$SRS.ID) ])
table(hmasm_samples$HMASM.QC)
#row.names(hmasm_samples) <- as.character(hmasm_samples$SRS.ID)

sites1 <- as.character(hmasm_samples$Body.Site)
names(sites1) <- as.character(hmasm_samples$SRS.ID)
sites2 <- as.character(hmiwgs_samples$Body.Site)
names(sites2) <- as.character(hmiwgs_samples$SRS.ID)

head(sites1)
head(sites2)

shared <- intersect(names(sites1),names(sites2))

if(!all(sites1[ shared ] == sites2[ shared ])){
  stop("ERROR: Not all samples match on body site")
}else{
  cat("All samples match\n")
}

row.names(hmasm_samples) <- as.character(hmasm_samples$SRS.ID)
hmiwgs_samples$HMASM <- as.character(hmiwgs_samples$SRS.ID) %in% as.character(hmasm_samples$SRS.ID)
hmiwgs_samples$HMASM.QC <- hmasm_samples[ as.character(hmiwgs_samples$SRS.ID), "HMASM.QC" ]
hmiwgs_samples$Reads.File.Location <- NULL
hmiwgs_samples$Reads.File.MD5 <- NULL
hmiwgs_samples$Reads.File.Size <- NULL
head(hmiwgs_samples)

# table(hmasm_samples[ as.character(hmiwgs_samples$SRS.ID), "HMASM.QC" ], useNA = "always")
# hmasm_samples[ as.character(all_runs$secondary_sample_accession), "Body.Site"]

row.names(hmiwgs_samples) <- as.character(hmiwgs_samples$SRS.ID)
head(all_runs)
all_runs$Body.Site <- hmiwgs_samples[ as.character(all_runs$secondary_sample_accession), "Body.Site"]
table(all_runs$Body.Site, useNA = "always")

all_runs$HMASM <- hmiwgs_samples[ as.character(all_runs$secondary_sample_accession), "HMASM"]
table(all_runs$HMASM, useNA = "always")

all_runs$HMASM.QC <- hmiwgs_samples[ as.character(all_runs$secondary_sample_accession), "HMASM.QC"]
table(all_runs$HMASM.QC, useNA = "always")


