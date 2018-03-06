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
#' @importFrom reshape2 melt
calculate_prevalence <- function(Dat, thres = 1, group = NULL){
  if(class(Dat) != "Dataset")
    stop("ERROR: A Dataset object must be passed", call. = TRUE )
  
  if(is.null(group)){
    group <- rep("All", length.out = ncol(Dat.bin$Tab))
    group_n <- table(group)
    varnames <- c("Group", "Taxon")
  }else{
    group <- Dat$Map[ , group ]
    group_n <- table(group)
    varnames <- c("Taxon", "Group")
  }
  
  Dat.bin <- create_dataset(Tab = 1*(Dat$Tab >= thres), Map = Dat$Map, Tax = Dat$Tax)
  Dat.bin <- pool_samples.default(Tab = Dat.bin$Tab, groups = group, FUN = sum)
  
  # res <- data.frame(Taxon = row.names(Dat$Tab), Count = rowSums(Dat$Tab >= thres))
  # res$Proportion <- res$Count / ncol(Dat$Tab)
  
  Res <- reshape2::melt(Dat.bin$Tab, varnames = varnames, value.name = "Count")
  Res <- Res[ , c("Taxon", "Group", "Count") ]
  
  if(length(group_n) == 1)
    Res$Group <- "All"
  
  Res$Proportion <- as.vector(Res$Count / group_n[ as.character(Res$Group) ])
  
  return(Res)
}
###################################

# Mapping file is the same everywhere
Map <- read.table(file = "~/micropopgen/data/hmp_16S/HMMCP/pds.metadata.bz2",
                  sep = "\t", header = TRUE)
row.names(Map) <- paste(Map$nap_id, Map$dataset, sep = ".")
# head(Map)

files <- data.frame(Name = c("hmmcp.v13.hq.otu",
                             "hmmcp.v13.hq.phylotype"),
                    counts = c("~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.otu.counts.bz2",
                               "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.counts.bz2"),
                    taxonomy = c("~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.otu.lookup.bz2",
                                 "~/micropopgen/data/hmp_16S/HMMCP/hmp1.v13.hq.phylotype.lookup.bz2"),
                    stringsAsFactors = FALSE)
files

#' Formats input form HMMCP 16S files
#' 
#' Only to be used to process 16S files from MOTHUR \
#' pipeline of HMP (HMMCP). Count files must be pre-edited
#' to remove trailing tab
#' 
#' @return A Dataset object
#' 
#' @author Sur Herrera Paredes
#' 
format_input <- function(name, counts_file, taxonomy_file, Map = Map,
                         collapse_level = NULL){
  # name <- files$Name[1]
  # counts_file <- files$counts[1]
  # taxonomy_file <- files$taxonomy[1]
  
  # Count table
  Tab <- read.table(file = counts_file,
                    sep = "\t", row.names = 1, header = TRUE)
  Tab <- t(Tab)
  # dim(Tab)
  # Tab[1:5,1:5]
  
  # Taxonomy
  Tax <- read.table(file = taxonomy_file,
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  colnames(Tax) <- c("ID", "Taxonomy")
  Tax$Taxonomy <- phylotype2rdp(Tax$Taxonomy)
  Tax$ID <- paste("X", Tax$ID, sep = "")
  row.names(Tax) <- Tax$ID
  # head(Tax)
  
  # Create dataset
  to_remove <- c("positive_control.PPS", "positive_control.may1",
                 "positive_gd.PPS", "positive_mock.PPS", 
                 "water_blank.PPS", "water_blank.may1")
  Tab <- Tab[ , setdiff(colnames(Tab), to_remove) ]
  
  if(length(setdiff(colnames(Tab), row.names(Map))) > 0)
    stop("ERROR1")
  if(length( setdiff(row.names(Tab), row.names(Tax))) > 0)
    stop("ERROR2")
 
  Map <- Map[ colnames(Tab), ]
  Tax <- Tax[ row.names(Tab), ]
  Dat <- create_dataset(Tab = Tab, Map = Map, Tax = Tax)
  
  if(!is.null(collapse_level)){
    Dat <- collapse_by_taxonomy(Dat = Dat, Group = NULL,
                                level = collapse_level, sepchar = ";", FUN = sum)
  }
 
  return(Dat)
}


i <- 1
Dat <- format_input(name = files$Name[i], counts_file = files$counts[i], 
                    taxonomy_file = files$taxonomy[i], Map = Map, collapse_level = 7)
# Subset sites
Dat <- subset(Dat, body_site %in% c("Buccal mucosa", "Supragingival plaque", "Tongue dorsum", "Stool"),
              clean = TRUE, drop = TRUE)

prev <- calculate_prevalence(Dat = Dat, thres = 1, group = "body_site")
prev$Taxonomy <- Dat$Tax[ as.character(prev$Taxon), "Taxonomy"]
head(prev)

# Sort
prev <- prev[ order(prev$Group, prev$Proportion , decreasing = TRUE), ]
prev$Taxon <- factor(prev$Taxon, unique(prev$Taxon))

# Plot
p1 <- ggplot(prev, aes(x = Taxon, y = Proportion, group = Group, col = Group )) +
  facet_wrap(~ Group, ncol = 1) +
  geom_line() +
  theme(axis.text.x = element_blank()) +
  theme_blackbox
p1
filename <- paste(files$Name[i], ".prevalence_by_site.svg", sep = "")
ggsave(filename, p1, width = 4, height = 8)



filename <- paste(files$Name[i], ".topprev.txt", sep = "")
write.table(subset(prev, Proportion >= 0.75), file = filename,
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)




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





