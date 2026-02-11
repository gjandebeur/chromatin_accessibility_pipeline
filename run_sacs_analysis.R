#!/usr/bin/env Rscript

################################################################################
# MICROEXON SACS ANALYSIS - COMPLETE PIPELINE
#
# This script runs the entire analysis from start to finish:
#   1. Data preparation & window creation
#   2. Chromatin mark processing (H3K27ac, H3K27me3, CTCF)
#   3. DNase processing (signal & peaks)
#   4. SACS score calculation
#   5. Differential splicing analysis
#   6. Figure generation
#   7. HOMER prep (BED file export)
#
# USAGE:
#   1. Edit config.txt with your file paths
#   2. Run: Rscript run_sacs_analysis.R
#
# Author: Gabe J, Wren Lab
################################################################################

cat("=== MICROEXON SACS ANALYSIS PIPELINE ===\n\n")

################################################################################
# LOAD CONFIG
################################################################################

cat("Loading configuration...\n")

# Read config file
config_lines <- readLines("config.txt")
config_lines <- config_lines[!grepl("^#", config_lines) & nchar(trimws(config_lines)) > 0]

config <- list()
for (line in config_lines) {
  parts <- strsplit(line, "=")[[1]]
  if (length(parts) == 2) {
    key <- trimws(parts[1])
    value <- trimws(gsub('"', '', parts[2]))
    
    # Try to convert to numeric if possible
    if (!is.na(suppressWarnings(as.numeric(value)))) {
      config[[key]] <- as.numeric(value)
    } else {
      config[[key]] <- value
    }
  }
}

cat(sprintf("  Input directory: %s\n", dirname(config$ME_ANNOTATIONS)))
cat(sprintf("  Output directory: %s\n", config$OUTPUT_DIR))

################################################################################
# LOAD LIBRARIES
################################################################################

cat("\nLoading libraries...\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(stringr)
})

# Create output directory
dir.create(config$OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
# STAGE 1: DATA PREPARATION
################################################################################

cat("\n=== STAGE 1: Data Preparation ===\n")

# Load microexon annotations
cat("Loading microexon annotations...\n")
ME_annot <- read_tsv(config$ME_ANNOTATIONS, show_col_types = FALSE) %>%
  separate(ME, into = c("chrom", "strand", "start", "end"),
           sep = "_", convert = TRUE, remove = FALSE) %>%
  mutate(ME_uid = paste(chrom, start, end, strand, sep = "_"),
         ME_length = end - start)

cat(sprintf("  Loaded %d microexons\n", nrow(ME_annot)))

# Load PSI values
cat("Loading PSI values and filtering...\n")
psi <- read_tsv(config$PSI_VALUES, show_col_types = FALSE) %>%
  mutate(original_coords = ME_coords) %>%
  separate_rows(ME_coords, sep = ",") %>%
  extract(ME_coords, into = c("chrom", "strand", "start", "end"),
          regex = "(chr[^_]+)_([+-])_([0-9]+)_([0-9]+)", convert = TRUE) %>%
  mutate(ME_uid = paste(chrom, start, end, strand, sep = "_")) %>%
  filter(!is.na(PSI), PSI > config$PSI_THRESHOLD)

cat(sprintf("  Filtered to %d high-confidence microexons (PSI > %.2f)\n",
            nrow(psi), config$PSI_THRESHOLD))

# Merge with annotations
ME_merged <- psi %>%
  left_join(ME_annot %>% select(ME_uid, ME_seq, ME_P_value, Transcript, ME_length),
            by = "ME_uid") %>%
  filter(!is.na(ME_seq))

# Create GRanges and windows
cat("Creating genomic windows...\n")
GR_ME <- GRanges(
  seqnames = ME_merged$chrom,
  ranges = IRanges(start = ME_merged$start, end = ME_merged$end),
  strand = ME_merged$strand,
  ME_uid = ME_merged$ME_uid,
  ME_seq = ME_merged$ME_seq,
  PSI = ME_merged$PSI
)

# Flanking windows
splice5 <- flank(GR_ME, width = config$FLANK_SIZE, start = TRUE)
splice3 <- flank(GR_ME, width = config$FLANK_SIZE, start = FALSE)
GR_ME$Window <- "ME_core"
splice5$Window <- "splice5"
splice3$Window <- "splice3"

all_windows <- c(GR_ME, splice5, splice3)

cat(sprintf("  Created %d windows (%d microexons × 3)\n",
            length(all_windows), length(GR_ME)))

# Save windows for reference
ME_windows_df <- as.data.frame(all_windows) %>%
  filter(!is.na(ME_seq))

write_tsv(ME_windows_df, file.path(config$OUTPUT_DIR, "ME_windows.tsv"))

################################################################################
# STAGE 2: CHROMATIN MARK PROCESSING
################################################################################

cat("\n=== STAGE 2: Chromatin Mark Processing ===\n")

# Helper function: extract BigWig signal
extract_bigwig_signal <- function(bw_dir, gr_object, name) {
  bw_files <- list.files(bw_dir, pattern = "\\.bw$|\\.bigwig$", full.names = TRUE)
  if (length(bw_files) == 0) {
    cat(sprintf("  WARNING: No BigWig files found in %s\n", bw_dir))
    return(NULL)
  }
  
  cat(sprintf("  Processing %d %s BigWig files...\n", length(bw_files), name))
  
  all_signals <- list()
  for (bw_file in bw_files) {
    bw <- import(bw_file, format = "BigWig")
    overlaps <- findOverlaps(gr_object, bw)
    
    if (length(overlaps) > 0) {
      signal_df <- data.frame(
        ME_uid = mcols(gr_object)$ME_uid[queryHits(overlaps)],
        signal = mcols(bw)$score[subjectHits(overlaps)]
      ) %>%
        group_by(ME_uid) %>%
        summarise(signal = mean(signal, na.rm = TRUE), .groups = 'drop')
      
      all_signals[[basename(bw_file)]] <- signal_df
    }
  }
  
  if (length(all_signals) == 0) return(NULL)
  
  bind_rows(all_signals) %>%
    group_by(ME_uid) %>%
    summarise(mean_signal = mean(signal, na.rm = TRUE), .groups = 'drop')
}

# Helper function: extract peak overlaps
extract_peak_overlaps <- function(bed_dir, gr_object, name) {
  bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
  if (length(bed_files) == 0) {
    cat(sprintf("  WARNING: No BED files found in %s\n", bed_dir))
    return(NULL)
  }
  
  cat(sprintf("  Processing %d %s BED files...\n", length(bed_files), name))
  
  all_overlaps <- list()
  for (bed_file in bed_files) {
    peaks <- read.table(bed_file, header = FALSE, fill = TRUE)
    if (ncol(peaks) < 3) next
    
    peaks_gr <- GRanges(seqnames = peaks[,1],
                       ranges = IRanges(start = peaks[,2], end = peaks[,3]))
    
    overlaps <- findOverlaps(gr_object, peaks_gr)
    has_peak <- as.integer(seq_along(gr_object) %in% queryHits(overlaps))
    
    all_overlaps[[basename(bed_file)]] <- data.frame(
      ME_uid = mcols(gr_object)$ME_uid,
      has_peak = has_peak
    )
  }
  
  if (length(all_overlaps) == 0) return(NULL)
  
  bind_rows(all_overlaps) %>%
    group_by(ME_uid) %>%
    summarise(has_peak = max(has_peak), .groups = 'drop')
}

# Process each chromatin mark
cat("\nProcessing H3K27ac (active mark)...\n")
H3K27ac_signal <- extract_bigwig_signal(config$H3K27AC_BIGWIG_DIR, GR_ME, "H3K27ac")
if (!is.null(H3K27ac_signal)) {
  threshold <- quantile(H3K27ac_signal$mean_signal, config$H3K27AC_TOP_PERCENTILE, na.rm = TRUE)
  H3K27ac_signal$SACS_H3K27ac <- as.integer(H3K27ac_signal$mean_signal >= threshold)
  cat(sprintf("  %d MEs in top 20%% (%.1f%%)\n",
              sum(H3K27ac_signal$SACS_H3K27ac),
              100 * mean(H3K27ac_signal$SACS_H3K27ac)))
  write_csv(H3K27ac_signal %>% select(ME_uid, SACS_H3K27ac),
            file.path(config$OUTPUT_DIR, "SACS_H3K27ac.csv"))
}

cat("\nProcessing H3K27me3 (repressive mark)...\n")
H3K27me3_peaks <- extract_peak_overlaps(config$H3K27ME3_BED_DIR, GR_ME, "H3K27me3")
if (!is.null(H3K27me3_peaks)) {
  H3K27me3_peaks$SACS_H3K27me3 <- -1 * H3K27me3_peaks$has_peak  # Negative for repressive
  cat(sprintf("  %d MEs with repressive peaks (%.1f%%)\n",
              sum(H3K27me3_peaks$has_peak),
              100 * mean(H3K27me3_peaks$has_peak)))
  write_csv(H3K27me3_peaks %>% select(ME_uid, SACS_H3K27me3),
            file.path(config$OUTPUT_DIR, "SACS_H3K27me3.csv"))
}

cat("\nProcessing CTCF (insulator)...\n")
CTCF_peaks <- extract_peak_overlaps(config$CTCF_BED_DIR, GR_ME, "CTCF")
if (!is.null(CTCF_peaks)) {
  CTCF_peaks$SACS_CTCF <- CTCF_peaks$has_peak
  cat(sprintf("  %d MEs with CTCF binding (%.1f%%)\n",
              sum(CTCF_peaks$SACS_CTCF),
              100 * mean(CTCF_peaks$SACS_CTCF)))
  write_csv(CTCF_peaks %>% select(ME_uid, SACS_CTCF),
            file.path(config$OUTPUT_DIR, "SACS_CTCF.csv"))
}

################################################################################
# STAGE 3: DNASE PROCESSING
################################################################################

cat("\n=== STAGE 3: DNase Processing ===\n")

cat("\nProcessing DNase signal...\n")
DNase_signal <- extract_bigwig_signal(config$DNASE_BIGWIG_DIR, GR_ME, "DNase")
if (!is.null(DNase_signal)) {
  threshold <- quantile(DNase_signal$mean_signal, config$DNASE_TOP_PERCENTILE, na.rm = TRUE)
  DNase_signal$SACS_DNase_signal <- as.integer(DNase_signal$mean_signal >= threshold)
  cat(sprintf("  %d MEs in top 20%% (%.1f%%)\n",
              sum(DNase_signal$SACS_DNase_signal),
              100 * mean(DNase_signal$SACS_DNase_signal)))
  write_csv(DNase_signal %>% select(ME_uid, mean_signal, SACS_DNase_signal),
            file.path(config$OUTPUT_DIR, "SACS_DNase_signal.csv"))
}

cat("\nProcessing DNase peaks...\n")
DNase_peaks <- extract_peak_overlaps(config$DNASE_BED_DIR, GR_ME, "DNase")
if (!is.null(DNase_peaks)) {
  DNase_peaks$SACS_DNase_peak <- DNase_peaks$has_peak
  cat(sprintf("  %d MEs with DNase peaks (%.1f%%)\n",
              sum(DNase_peaks$SACS_DNase_peak),
              100 * mean(DNase_peaks$SACS_DNase_peak)))
  write_csv(DNase_peaks %>% select(ME_uid, SACS_DNase_peak),
            file.path(config$OUTPUT_DIR, "SACS_DNase_peak.csv"))
}

################################################################################
# STAGE 4: SACS INTEGRATION
################################################################################

cat("\n=== STAGE 4: SACS Integration ===\n")

# Start with all MEs
SACS_all <- data.frame(ME_uid = mcols(GR_ME)$ME_uid)

# Merge all components
if (!is.null(DNase_signal)) {
  SACS_all <- SACS_all %>%
    left_join(DNase_signal %>% select(ME_uid, SACS_DNase_signal), by = "ME_uid")
}
if (!is.null(DNase_peaks)) {
  SACS_all <- SACS_all %>%
    left_join(DNase_peaks %>% select(ME_uid, SACS_DNase_peak), by = "ME_uid")
}
if (!is.null(H3K27ac_signal)) {
  SACS_all <- SACS_all %>%
    left_join(H3K27ac_signal %>% select(ME_uid, SACS_H3K27ac), by = "ME_uid")
}
if (!is.null(H3K27me3_peaks)) {
  SACS_all <- SACS_all %>%
    left_join(H3K27me3_peaks %>% select(ME_uid, SACS_H3K27me3), by = "ME_uid")
}
if (!is.null(CTCF_peaks)) {
  SACS_all <- SACS_all %>%
    left_join(CTCF_peaks %>% select(ME_uid, SACS_CTCF), by = "ME_uid")
}

# Replace NA with 0 and calculate total
SACS_all <- SACS_all %>%
  mutate(
    across(starts_with("SACS_"), ~replace_na(.x, 0)),
    SACS_total = rowSums(select(., starts_with("SACS_")), na.rm = TRUE)
  )

# Add annotations
SACS_all <- SACS_all %>%
  left_join(ME_annot %>% select(ME_uid, chrom, start, end, strand, ME_length, Transcript),
            by = "ME_uid") %>%
  left_join(ME_merged %>% select(ME_uid, PSI, ME_seq), by = "ME_uid")

write_csv(SACS_all, file.path(config$OUTPUT_DIR, "SACS_all_annotated.csv"))

cat("\nSACS Score Distribution:\n")
print(table(SACS_all$SACS_total))

cat(sprintf("\nHigh-scoring MEs (SACS = 4): %d (%.1f%%)\n",
            sum(SACS_all$SACS_total == 4, na.rm = TRUE),
            100 * mean(SACS_all$SACS_total == 4, na.rm = TRUE)))

################################################################################
# STAGE 5: DIFFERENTIAL SPLICING (OPTIONAL)
################################################################################

if (file.exists(config$DIFF_SPLICING)) {
  cat("\n=== STAGE 5: Differential Splicing Analysis ===\n")
  
  dpsi <- read_tsv(config$DIFF_SPLICING, show_col_types = FALSE) %>%
    mutate(direction = case_when(
      DeltaPsi >= config$DPSI_THRESHOLD ~ "AD_up",
      DeltaPsi <= -config$DPSI_THRESHOLD ~ "AD_down",
      TRUE ~ "NS"
    )) %>%
    filter(direction != "NS")
  
  write_csv(dpsi, file.path(config$OUTPUT_DIR, "differential_splicing_filtered.csv"))
  cat(sprintf("  %d MEs with |dPSI| >= %.2f\n", nrow(dpsi), config$DPSI_THRESHOLD))
}

################################################################################
# STAGE 6: GENERATE FIGURES
################################################################################

cat("\n=== STAGE 6: Generating Figures ===\n")

# SACS distribution
p1 <- ggplot(SACS_all, aes(x = factor(SACS_total))) +
  geom_bar(fill = "steelblue", alpha = 0.8, color = "black", width = 0.7) +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5) +
  labs(title = "SACS Score Distribution",
       subtitle = sprintf("n = %d microexons", nrow(SACS_all)),
       x = "SACS Score", 
       y = "Number of Microexons") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))

ggsave(file.path(config$OUTPUT_DIR, "SACS_distribution.pdf"),
       p1, width = 8, height = 6)
cat("  Saved: SACS_distribution.pdf\n")

# Component summary
component_counts <- SACS_all %>%
  summarise(
    `DNase Signal` = sum(SACS_DNase_signal, na.rm = TRUE),
    `DNase Peak` = sum(SACS_DNase_peak, na.rm = TRUE),
    `H3K27ac` = sum(SACS_H3K27ac, na.rm = TRUE),
    `H3K27me3\n(repressive)` = sum(SACS_H3K27me3 == -1, na.rm = TRUE),
    `CTCF` = sum(SACS_CTCF, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "Component", values_to = "Count")

p2 <- ggplot(component_counts, aes(x = reorder(Component, Count), y = Count)) +
  geom_col(fill = "coral", alpha = 0.8, color = "black", width = 0.7) +
  geom_text(aes(label = Count), hjust = -0.2, size = 4) +
  coord_flip() +
  labs(title = "SACS Component Distribution",
       subtitle = "Number of microexons with each chromatin feature",
       x = "", 
       y = "Number of Microexons") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12))

ggsave(file.path(config$OUTPUT_DIR, "SACS_components.pdf"),
       p2, width = 8, height = 6)
cat("  Saved: SACS_components.pdf\n")

################################################################################
# STAGE 7: PREPARE FOR HOMER
################################################################################

cat("\n=== STAGE 7: Preparing BED files for HOMER ===\n")

# Extract high-SACS microexons
high_SACS <- SACS_all %>%
  filter(SACS_total >= 3, !is.na(chrom)) %>%
  head(config$HOMER_TOP_N) %>%
  mutate(
    bed_start = start - config$HOMER_WINDOW_SIZE,
    bed_end = end + config$HOMER_WINDOW_SIZE,
    score = 0
  ) %>%
  select(chrom, bed_start, bed_end, ME_uid, score, strand)

write.table(high_SACS,
            file.path(config$OUTPUT_DIR, "high_SACS_for_HOMER.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("  Created BED file with %d high-SACS regions\n", nrow(high_SACS)))
cat(sprintf("  Ready for HOMER: findMotifsGenome.pl high_SACS_for_HOMER.bed %s output/\n",
            config$HOMER_GENOME))

################################################################################
# COMPLETE
################################################################################

cat("\n=== PIPELINE COMPLETE ===\n\n")
cat("Results saved to:", config$OUTPUT_DIR, "\n\n")
cat("Output files:\n")
cat("  - SACS_all_annotated.csv         (main results)\n")
cat("  - SACS_distribution.pdf          (score distribution)\n")
cat("  - SACS_components.pdf            (component breakdown)\n")
cat("  - high_SACS_for_HOMER.bed        (for motif analysis)\n")
if (file.exists(config$DIFF_SPLICING)) {
  cat("  - differential_splicing_filtered.csv\n")
}
cat("\nNext steps:\n")
cat("  1. Review SACS_all_annotated.csv for high-scoring microexons\n")
cat("  2. Run HOMER on high_SACS_for_HOMER.bed\n")
cat("  3. Validate results with your controls\n\n")
