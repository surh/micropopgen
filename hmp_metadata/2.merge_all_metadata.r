# Load metadata from metagenome samples
# file <- "~/micropopgen/data/HMP/HMIWGS_healthy_run_metadata/all_runs.txt"
# all_runs <- read.table(file,header = TRUE,sep = "\t")
# head(all_runs)

# file <- "~/micropopgen/data/HMP/HMIWGS_healthy.csv"
# hmiwgs_samples <- read.table(file,sep = ",", quote = '"', header = TRUE)
# head(hmiwgs_samples)

file <- "~/micropopgen/data/HMP/HMASM.csv"
hmasm_samples <- read.table(file,sep = ",", quote = '"', header = TRUE)
head(hmasm_samples)

file <- "~/micropopgen/data/HMP/HMASM-690.csv"
hmasm_samples_690 <- read.table(file,header = TRUE,sep = ",",quote = '"')
head(hmasm_samples_690)

file <- "~/micropopgen/data/HMP/hmp_catalogue_wgs_all_metadata.csv"
catalogue <- read.csv(file)
head(catalogue)

# file <- "subject_metadata.txt"
# ssu <- read.table(file,header = TRUE,sep = "\t",quote = '')
# head(ssu)

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


# Add body subsite
hmiwgs_samples$Body.Subsite <- NA
hmiwgs_samples$Body.Subsite[ hmiwgs_samples$Body.Site %in% c("anterior_nares") ] <- "Airways"
hmiwgs_samples$Body.Subsite[ hmiwgs_samples$Body.Site %in% c("stool") ] <- "Gastrointestinal_tract"
hmiwgs_samples$Body.Subsite[ hmiwgs_samples$Body.Site %in% c("attached_keratinized_gingiva","buccal_mucosa","hard_palate",
                                                             "palatine_tonsils","saliva","aubgingival_plaque","supragingival_plaque",
                                                             "throat","tongue_dorsum","subgingival_plaque") ] <- "Oral"
hmiwgs_samples$Body.Subsite[ hmiwgs_samples$Body.Site %in% c("left_antecubital_fossa","left_retroauricular_crease","right_antecubital_fossa",
                                                             "right_retroauricular_crease") ] <- "Skin"
hmiwgs_samples$Body.Subsite[ hmiwgs_samples$Body.Site %in% c("mid_vagina","posterior_fornix","vaginal_introitus") ] <- "Urogenital_tract"
table(hmiwgs_samples$Body.Subsite,useNA = "always")
subset(hmiwgs_samples,is.na(Body.Subsite))

# Add metadata to run file
row.names(hmiwgs_samples) <- as.character(hmiwgs_samples$SRS.ID)
head(all_runs)
all_runs$Body.Site <- hmiwgs_samples[ as.character(all_runs$secondary_sample_accession), "Body.Site"]
table(all_runs$Body.Site, useNA = "always")

all_runs$Body.Subsite <- hmiwgs_samples[ as.character(all_runs$secondary_sample_accession), "Body.Subsite"]
table(all_runs$Body.Subsite, useNA = "always")

all_runs$HMASM <- hmiwgs_samples[ as.character(all_runs$secondary_sample_accession), "HMASM"]
table(all_runs$HMASM, useNA = "always")

all_runs$HMASM.QC <- hmiwgs_samples[ as.character(all_runs$secondary_sample_accession), "HMASM.QC"]
table(all_runs$HMASM.QC, useNA = "always")

# check no samples with infor are missing subject id
intersect(as.character((subset(all_runs, !is.na(subject_id))$secondary_sample_accession)),as.character((subset(all_runs, is.na(subject_id))$secondary_sample_accession)))

write.table(all_runs, "merged_metadata.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
