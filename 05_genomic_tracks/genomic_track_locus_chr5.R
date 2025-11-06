# Load required libraries
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
library(data.table)
library(dplyr)
library(ggplot2)         # Added for theme_set
library(extrafont)  # For managing fonts
library(AnnotationDbi)   # For mapIds
# Set working directory (adjust accordingly)
setwd("~/Projects/cfChIP")
source("../Scripts/functions.R")
theme_set(
  theme_classic(base_family = "Helvetica", base_size = 7)
)


# Define genomic region
loc <- "chr5:1889346-1889346"
win = 4e3


location <- get_loc(loc)
chromosome <- location$chromosome
start <- location$start -win
end <- location$end + win/2

# Define paths to BigWig files
wbk27.bw_file <- list.files('Results/bw_files/wbk27', full.names = TRUE)
k27.bw_file <- list.files('Results/bw_files/k27', full.names = TRUE, recursive = TRUE)

# Get global max for all the cfChIP and WB BigWig files
k27_max <- get_global_max(k27.bw_file, chromosome, start, end)
wb_max <- get_global_max(wbk27.bw_file, chromosome, start, end)
global_max <- max(k27_max, wb_max)

# Define improved color scheme
gwas_color <- "#D73027"   # Bright Red for GWAS significance
cfchip_color <- "blue" # Deep Blue for cfChIP-seq
h3k27ac_color <- "grey" # Light Blue for H3K27ac in whole blood
snp_color <- "#E31A1C"    # Dark Red for SNPs
gene_color <- "#1A9850"   # Deep Green for Gene annotation

# Create overlay tracks for each condition
k27_track <- create_data_track(k27.bw_file, chromosome, start, end, cfchip_color, "cfChIP", global_max)
wb_track <- create_data_track(wbk27.bw_file, chromosome, start, end, h3k27ac_color, "WBC ChIP", global_max)

# Gene annotation track
customFromTxDb <- GeneRegionTrack(
  TxDb.Hsapiens.UCSC.hg19.knownGene,
  chromosome = chromosome,
  geneSymbol = TRUE,
  transcriptAnnotation = "symbol",
  name = "Gene",
  background.panel = "white",
  background.title = "lightgray",
  col = gene_color,
  fill = adjustcolor(gene_color, alpha.f = 0.6)
)

# Map symbols using Homo.sapiens
z <- ranges(customFromTxDb)
z$symbol <- mapIds(Homo.sapiens, z$symbol, "SYMBOL", "TXNAME", multiVals = "first")
ranges(customFromTxDb) <- z

# GWAS significance data
gwas_data <- fread('Data/PRAD.GWAS.sumstats.txt')
gwas_data <- gwas_data %>%
  filter(chromosome == chromosome, base_pair_location <= end, base_pair_location > start)

gwas_track <- DataTrack(
  data = -log10(gwas_data$p_value),
  start = gwas_data$base_pair_location,
  end = gwas_data$base_pair_location,
  chromosome = chromosome,
  name = "GWAS\n-log10(P)",
  genome = "hg19",
  type = "p",
  ylim = c(0, max(-log10(gwas_data$p_value)) + 1),
  col = gwas_color,
  background.panel = "white",
  background.title = "lightgray",
  cex = 0.5
)

# MNLP annotation track
##Reference Sup data 5 of PMID: 37612286
MNLP_gr <- GRanges(
  seqnames = 5,             # Chromosome
  ranges = IRanges(start = 1889299, end = 1889346),  # SNP positions
  strand = "*"                         # SNPs are unstranded
)

snp_track <- DataTrack(
  start = start(MNLP_gr),
  end = start(MNLP_gr),
  chromosome = chromosome,
  genome = "hg19",
  name = "Indel",
  data = rep(1, length(MNLP_gr)),
  type = "h",
  col = snp_color,
  lwd = 1.5,
  ylim = c(0, 1),  # Adjust to match gene annotation track height
  showAxis=FALSE,
  background.panel = "white",
  background.title = "lightgray"
)

# Genome axis for reference
genome_axis <- GenomeAxisTrack(cex=0.5, 
                               background.title = "white")
itrack <- IdeogramTrack(genome = genome(customFromTxDb), chromosome = chromosome, cex=0.5)
# Combine all tracks
#all_tracks <- list(itrack, genome_axis, gwas_track, k27_track, wb_track, snp_track, customFromTxDb)
all_tracks <- list(genome_axis, gwas_track, k27_track, wb_track, snp_track, customFromTxDb)
all_tracks <- lapply(all_tracks, function(trk) {
  displayPars(trk) <- list(innerMargin = c(1.5, 0, 0, 0))  # bottom, left, top, right
  trk
})
# Plot the tracks
pdf(
  paste0("Manuscripts/Figures/1_f_raw.pdf"),
  width = 120/25.4, height = 100/25.4,   # mm -> inches
  family = "Helvetica", useDingbats = FALSE
)

plotTracks(
  all_tracks,
  from = start, to = end,
  sizes = c(3.5, 4, 4, 4, 2, 2),
  col.axis = "black", col.title = "black", fontcolor = "black",
  fontfamily = "Helvetica",
  cex = 10/12,           # ~10 pt body text
  cex.title = 10/12,     # ~10 pt titles
  cex.axis = 5.5/12,      # ~10 pt axes
  margin = 0,
  #innerMargin = c(2, 0, 0, 0),
  collapseTranscripts = "longest",
  background.panel = "white",
  title.width = 1,
  fontface.title = 1
)

dev.off()
