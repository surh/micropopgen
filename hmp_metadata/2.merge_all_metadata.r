# Load metadata from metagenome samples
file <- "~/micropopgen/data/HMP/2017-07-17.all_runs/all_runs.txt"
all_runs <- read.table(file,header = TRUE,sep = "\t")
head(all_runs)

file <- "~/micropopgen/data/HMP/HMASM.csv"
hmasm_samples <- read.table(file,sep = ",", quote = '"', header = TRUE)
head(hmasm_samples)

file <- "~/micropopgen/data/HMP/HMASM-690.csv"
hmasm_samples_690 <- read.table(file,header = TRUE,sep = ",",quote = '"')
head(hmasm_samples_690)

file <- "~/micropopgen/data/HMP/hmp_catalogue_wgs_all_metadata.csv"
catalogue <- read.csv(file)
head(catalogue)

# The first thing is to identify which samples are in the original publication,
# and which passed QC

# Compare body sites from hmasm sets
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
head(hmasm_samples)
#row.names(hmasm_samples) <- as.character(hmasm_samples$SRS.ID)

# Add information about puiblication and QC into main table
ftable(HMP.Isolation.Body.Site ~ HMP.Isolation.Body.Subsite, catalogue)
row.names(hmasm_samples) <- as.character(hmasm_samples$SRS.ID)
catalogue$HMASM <- as.character(catalogue$Sequence.Read.Archive.ID) %in% row.names(hmasm_samples)
catalogue$HMASM.QC <- hmasm_samples[ as.character(catalogue$Sequence.Read.Archive.ID), "HMASM.QC" ]
table(catalogue$HMASM, useNA = "always")
table(catalogue$HMASM.QC, useNA = "always")

all(as.character(all_runs$secondary_sample_accession) %in% as.character(catalogue$Sequence.Read.Archive.ID))

# Set to stool all GI samples
catalogue$HMP.Isolation.Body.Subsite <- as.character(catalogue$HMP.Isolation.Body.Subsite)
catalogue$HMP.Isolation.Body.Subsite[ catalogue$HMP.Isolation.Body.Site == "gastrointestinal_tract" ] <- "Stool"

# Now we need to identify samples for basic analysis
dat <- subset(catalogue, HMASM.QC == TRUE & HMP.Isolation.Body.Site %in% c("gastrointestinal_tract","oral"))
dat <- droplevels(dat)

dat.runs <- all_runs
# Get illumina runs only
#dat.runs <- subset(all_runs, instrument_model == "Illumina HiSeq 2000")

dat.runs <- droplevels(dat.runs)
table(dat.runs$instrument_model, useNA = "always")
table(dat.runs$library_strategy, useNA = "always")
table(dat.runs$library_source, useNA = "always")
table(dat.runs$library_layout, useNA = "always")
table(dat.runs$library_selection, useNA = "always")

# Get runs from selected samples
dat.runs <- subset(dat.runs, secondary_sample_accession %in% as.character(dat$Sequence.Read.Archive.ID))

dat.runs <- droplevels(dat.runs)
table(dat.runs$instrument_model, useNA = "always")
table(dat.runs$library_strategy, useNA = "always")
table(dat.runs$library_source, useNA = "always")
table(dat.runs$library_layout, useNA = "always")
table(dat.runs$library_selection, useNA = "always")

# Eliminate 454
dat.runs <- subset(dat.runs, instrument_model != "454 GS FLX Titanium")

dat.runs <- droplevels(dat.runs)
table(dat.runs$instrument_model, useNA = "always")
table(dat.runs$library_strategy, useNA = "always")
table(dat.runs$library_source, useNA = "always")
table(dat.runs$library_layout, useNA = "always")
table(dat.runs$library_selection, useNA = "always")

# Eliminate unpaired
dat.runs <- subset(dat.runs, library_layout == "PAIRED")

dat.runs <- droplevels(dat.runs)
table(dat.runs$instrument_model, useNA = "always")
table(dat.runs$library_strategy, useNA = "always")
table(dat.runs$library_source, useNA = "always")
table(dat.runs$library_layout, useNA = "always")
table(dat.runs$library_selection, useNA = "always")

# Select samples with runs
dat <- subset(dat, Sequence.Read.Archive.ID %in% unique(as.character(dat.runs$secondary_sample_accession)))

# Add depth
reads <- aggregate(read_count ~ secondary_sample_accession, data = dat.runs, FUN = sum)
summary(reads)
row.names(reads) <- as.character(reads$secondary_sample_accession)
dat$Reads <- reads[ as.character(dat$Sequence.Read.Archive.ID), "read_count"]

# Select samples with one visit and one replicate
dat <- subset(dat, Visit.Number == 1)
dat <- subset(dat, Replicate.Number == 0)

table(dat$HMP.Isolation.Body.Site, useNA = "always")
table(dat$HMP.Isolation.Body.Subsite, useNA = "always")

# Select samples with subjects that donated multiple samples
subjects <- table(dat$MRN..Subject.ID)
subjects <- names(subjects)[subjects > 1]
dat <- subset(dat, dat$MRN..Subject.ID %in% subjects)
dat <- droplevels(dat)

table(dat$HMP.Isolation.Body.Site, useNA = "always")
table(dat$HMP.Isolation.Body.Subsite, useNA = "always")

# Select samples in 4 more common sites
dat <- subset(dat,HMP.Isolation.Body.Subsite %in% c("Stool","Supragingival plaque","Buccal mucosa","Tongue dorsum"))
subjects <- table(dat$MRN..Subject.ID)
subjects <- names(subjects)[subjects == 4 ]
dat <- subset(dat, dat$MRN..Subject.ID %in% subjects)
dat <- droplevels(dat)

table(dat$HMP.Isolation.Body.Site, useNA = "always")
table(dat$HMP.Isolation.Body.Subsite, useNA = "always")

ftable(dat$Host.Gender ~ dat$HMP.Isolation.Body.Subsite)

write.table(dat,"samples_to_start.txt",sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## Select runs
dat.runs <- subset(dat.runs, secondary_sample_accession %in% as.character(dat$Sequence.Read.Archive.ID))
dat.runs <- droplevels(dat.runs)
download <- data.frame(Sample = dat.runs$secondary_sample_accession, Run = dat.runs$run_accession)

download <- download[ order(download$Sample), ]
head(download)

write.table(dat,"runs_to_download.txt",sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


