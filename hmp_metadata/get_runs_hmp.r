library(SRAdb)

# Set up connection
sqlfile <-'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
dbcon <- dbConnect(SQLite(),sqlfile)

# sraConvert(c('SRS049712','SRS049995'),'run',dbcon)

# Read tables
hmasm <- read.csv("~/micropopgen/data/HMP/HMASM.csv")
hmasm.690 <- read.csv("~/micropopgen/data/HMP/HMASM-690.csv")
hmiwgs <- read.csv("~/micropopgen/data/HMP/HMIWGS_healthy.csv")

# Check consistency
all(hmasm.690$SRS.ID %in% hmasm$SRS.ID)
all(hmasm$SRS.ID %in% hmiwgs$SRS.ID)

# Get runs
res <- sraConvert(hmiwgs$SRS.ID, c("sra"),dbcon)

write.table(res, "hmiwgs_run_table.txt", quote = FALSE, col.names = TRUE, row.names = FALSE,sep = "\t")
