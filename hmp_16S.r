library(AMOR)

########## FUNCTIONS #################
phylotype2rdp <- function(x, split.char = ';'){
  for(i in 1:length(x)){
    x[i] <- gsub(pattern = "\\(\\d+\\)", replacement = "", x = x[i], perl = TRUE)
  }
  
  return(x)
}

phylotype2rdp2 <- function(x, split.char = ';'){
  sapply(strsplit(x, split = split.char), function(x){
    
    # Get last quality
    # last <- x[ length(x) ]
    # pat <- "\\((\\d+)\\)"
    # match <- regexpr(pattern = pat, text = last, perl = TRUE)
    # start.match <- attr(x = match, which = "capture.start")
    # length.match <- attr(x = match, which = "capture.length")
    # q <- substr(x = last, start = start.match, stop = start.match + length.match - 1)
    #return(q)
    
    # Clean
    x <- paste(gsub(pattern = pat, replacement = "", x = x), collapse = ";")
    return(x)
  })
  
}
###################################
Tab <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.counts.bz2",
                  sep = "\t", row.names = 1, header = TRUE)
dim(Tab)
Tab[1:5,1:5]

Map <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/pds.metadata.bz2",
                   sep = "\t", header = TRUE)
row.names(Map) <- paste(Map$nap_id, Map$dataset, sep = ".")
head(Map)

Tax <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.lookup.bz2",
                  sep = "\t", header = TRUE)
head(Tax)

# Proces taxonomy file
colnames(Tax) <- c("ID","Taxonomy")




# x <- as.character(Tax$Taxonomy[1:200])
# x
# system.time(phylotype2rdp(x = x))
# system.time(phylotype2rdp2(x = x))
# 
# dat <- NULL
# for(f in seq(from = 100, to = 600, by = 100)){
#   x <- as.character(Tax$Taxonomy[1:f])
#   
#   
#   t1 <- system.time(phylotype2rdp(x = x))
#   t2 <- system.time(phylotype2rdp2(x = x))
#   
#   t1 <- rbind(t1)
#   colnames(t1) <- paste("t1", colnames(t1), sep = ".")
#   row.names(t1) <- as.character(f)
#   
#   t2 <- rbind(t2)
#   colnames(t2) <- paste("t2", colnames(t2), sep = ".")
#   row.names(t2) <- as.character(f)
#   
#   dat <- rbind(dat,cbind(t1, t2))
# }
# dat





Dat <- create_dataset(Tab = t(Tab))