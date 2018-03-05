library(AMOR)

########## FUNCTIONS #################
#' Format HMMCP phylotype taxonomy for AMIR
#' 
#' @param x A character vector, one element per taxa
#' @param split.char A character string indicating the field
#' delimiter character
#' 
#' @return A character vector
#' 
#' @author Sur Herrera Paredes
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
# Count table
Tab <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.counts.bz2",
                  sep = "\t", row.names = 1, header = TRUE)
Tab <- t(Tab)
dim(Tab)
Tab[1:5,1:5]

# Mapping file
Map <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/pds.metadata.bz2",
                   sep = "\t", header = TRUE)
row.names(Map) <- paste(Map$nap_id, Map$dataset, sep = ".")
head(Map)

# Taxonomy
Tax <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.lookup.bz2",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(Tax) <- c("ID", "Taxonomy")
Tax$Taxonomy <- phylotype2rdp(Tax$Taxonomy)
Tax$ID <- paste("X", Tax$ID, sep = "")
row.names(Tax) <- Tax$ID
head(Tax)

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

# Create dataset
to_remove <- c("positive_control.PPS", "positive_control.may1",
               "positive_gd.PPS", "positive_mock.PPS", 
             "water_blank.PPS", "water_blank.may1")
Tab <- Tab[ , setdiff(colnames(Tab), to_remove) ]

setdiff(colnames(Tab), row.names(Map))
setdiff(row.names(Tab), row.names(Tax))

Map <- Map[ colnames(Tab), ]
Tax <- Tax[ row.names(Tab), ]
Dat <- create_dataset(Tab = Tab, Map = Map, Tax = Tax)
Dat

# Subset sites
Dat <- subset(Dat, body_site %in% c("Buccal mucosa", "Supragingival plaque", "Tongue dorsum", "Stool"),
              clean = TRUE, drop = TRUE)


#' Calculate taxon prevalence
#' 
#' Calculates the prevalence of each
#' taxon overall or by some grouping factor.
#' 
#' @param Dat a dataset object
#' @param thres Minimum number of reads for a taxon
#' in a sample to be counted as present.
#' @param group A grouping variable
#' 
#' @author Sur Herrera Paredes
#' 
calculate_prevalence <- function(Dat, thres = 1, group = NULL){
  if(class(Dat) != "Dataset")
    stop("ERROR: A Dataset object must be passed", call. = TRUE )
  
}


