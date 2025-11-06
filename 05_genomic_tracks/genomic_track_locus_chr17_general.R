# 🔹 Load Required Libraries
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)

#Set Working Directory (Adjust Accordingly)
setwd("~/Projects/cfChIP")

#Define Gene Name and Genomic Location
Gene.name <- 'HOXB13'  # Update this to the relevant gene
loc <- "chr17:46802126-46806112"

#Define Paths to BigWig Files
wbk27.bw_file <- list.files('Results/bw_files/wbk27', full.names = TRUE)
k27.NEPC.bw_file <- list.files('Results/bw_files/k27/NEPC/', full.names = TRUE, recursive = TRUE)
k27.CRC.bw_file <- list.files('Results/bw_files/k27/colorectal/', full.names = TRUE, recursive = TRUE)
k27.BRCA.bw_file <- list.files('Results/bw_files/k27/Breast/', full.names = TRUE, recursive = TRUE)
k27.PRAD.bw_file <- list.files('Results/bw_files/k27/PRAD/', full.names = TRUE, recursive = TRUE)

#Function to Extract Chromosome, Start, and End Positions
get_loc <- function(loc_string) {
  match <- regmatches(loc_string, regexpr("chr[0-9XYM]+", loc_string))
  chromosome <- match
  positions <- strsplit(gsub("chr[0-9XYM]+:", "", loc_string), "-")[[1]]
  start <- as.numeric(positions[1])
  end <- as.numeric(positions[2])
  return(list(chromosome = chromosome, start = start, end = end))
}

#Extract Location Data
location <- get_loc(loc)
chromosome <- location$chromosome
start <- location$start
end <- location$end

#Function to Compute Global Maximum for Consistent ylim
get_global_max <- function(bw_files, chromosome, start, end) {
  max_values <- sapply(bw_files, function(bw) {
    bw_data <- import(bw, which = GRanges(chromosome, IRanges(start, end)))
    max(score(bw_data), na.rm = TRUE)
  })
  return(max(max_values, na.rm = TRUE))
}

#Compute Global Max for Consistency Across Tracks
NEPC_max <- get_global_max(k27.NEPC.bw_file, chromosome, start, end)
BRCA_max <- get_global_max(k27.BRCA.bw_file, chromosome, start, end)
PRAD_max <- get_global_max(k27.PRAD.bw_file, chromosome, start, end)
CRC_max <- get_global_max(k27.CRC.bw_file, chromosome, start, end)
wb_max <- get_global_max(wbk27.bw_file, chromosome, start, end)

global_max <- max(NEPC_max, PRAD_max, BRCA_max, CRC_max, wb_max)

#Function to Create DataTrack with Consistent ylim
create_data_track <- function(bw_files, chromosome, start, end, color, track_name, global_max) {
  # Create a DataTrack for each BigWig file and combine them
  data_tracks <- lapply(bw_files, function(bw_file) {
    short_name <- gsub(".rep1_treat_pileup.bw", "", basename(bw_file))
    short_name <- substr(short_name, 1, 15)  # Limit to 15 characters
    
    DataTrack(range = bw_file, 
              genome = "hg19", 
              chromosome = chromosome, 
              name = track_name, 
              type = "h", 
              ylim = c(0, global_max + global_max * 0.1),  # Set consistent ylim
              col.histogram = color, 
              fill.histogram = color, 
              lwd = 2, 
              cex.title = 0.8)
  })
  # Create an OverlayTrack to combine them
  OverlayTrack(trackList = data_tracks, name = track_name)
}
# 🔹 Create DataTracks for Each Condition
k27.BRCA_track <- create_data_track(k27.BRCA.bw_file, chromosome, start, end, 'grey', "BRCA", global_max)
k27.CRC_track <- create_data_track(k27.CRC.bw_file, chromosome, start, end, 'grey', "CRC", global_max)
k27.NEPC_track <- create_data_track(k27.NEPC.bw_file, chromosome, start, end, 'grey', "NEPC", global_max)
k27.PRAD_track <- create_data_track(k27.PRAD.bw_file, chromosome, start, end, 'blue', "PRAD", global_max)
wb_track <- create_data_track(wbk27.bw_file, chromosome, start, end, 'grey', "WB H3K27ac", global_max)

# 📌 Gene Annotation Track
customFromTxDb <- GeneRegionTrack(
  TxDb.Hsapiens.UCSC.hg19.knownGene, 
  chromosome = chromosome, 
  geneSymbol = TRUE, 
  transcriptAnnotation = "symbol", 
  name = "Gene", 
  background.panel = "white", 
  background.title = "grey",
  col = "#1A9850",
  fill = adjustcolor("#1A9850", alpha.f = 0.6)
)

#Map Symbols Using Homo.sapiens
z <- ranges(customFromTxDb)
z$symbol <- mapIds(Homo.sapiens, z$symbol, "SYMBOL", "TXNAME", multiVals = "first")
ranges(customFromTxDb) <- z

#Genome Axis Track
genome_axis <- GenomeAxisTrack()

#Combine Tracks (Excluding Gene Track if Not Needed)
all_tracks.nogene <- list(genome_axis, k27.CRC_track, k27.NEPC_track, k27.PRAD_track, wb_track, customFromTxDb)

#Save and Plot the Main Figure
png(paste0("Figures/Gviz/", Gene.name, '.png'), width = 16, height = 32, res = 200, units = "cm")
plotTracks(
  all_tracks.nogene, 
  from = start, to = start+2000,  
  main = "Overlay of cfChIP and WB ChIP H3K27ac", 
  cex.title = 0.7, 
  fontsize = 10, 
  add53 = TRUE,  
  title.width = 0.6, 
  type = "hist"
)
dev.off()

#Save and Plot the Gene Label Track
png(paste0("Figures/Gviz/", Gene.name, '.label.png'), width = 16, height = 8, res = 200, units = "cm")
plotTracks(
  customFromTxDb, 
  from = start, to = start+2000,   
  main = "Overlay of cfChIP and WB ChIP H3K27ac", 
  cex.title = 0.7, 
  fontsize = 10, 
  add53 = TRUE,  
  title.width = 0.6, 
  background.panel = "white"
  )
dev.off()
