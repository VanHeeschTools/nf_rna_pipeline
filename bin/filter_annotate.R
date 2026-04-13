#!/usr/bin/env Rscript

# Input requirements: =========================================================
# 1. gtf_ref_loc: GTF file with ensembl annotations, currently using hg38 v 102
# 2. gtf_novel_file: generated using gffcompare to merge all stringtie output GTF files
# 3. tracking_file: generated in same way as gtf_novel_file, contains mappings between XLOC / TCONS IDs and reference ENSG and ENST IDs
# 4. min_occurrence: minimum number of samples in which a transcript must be found to be included in the final GTF
# 5. min_tpm: minimum expression value in which a transcript must be found to be included in the final GTF
# 6. outfile: output prefix
# 7. scripts_dir: directory containing filter_annotate_functions.R
# 8. (OPTIONAL) gtf_refseq_basename: GTF with refseq annotations, used to annotate XR and NR overlap

#
# Output: =====================================================================
# 1. outfile.gtf: filtered custom annotated GTF
# 2. outfile.log: logging file with version info and filtering steps
# 1. outfile.tsv: expression matrix of novel transcripts (from StringTie)

# Track versions:
script_version <- "2.1"
script_date <- "08-12-2025"

# Load required packages
message(paste(Sys.time(), "Loading required libraries ..."), sep = "\t")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# Global arguments
args <- commandArgs(trailingOnly = TRUE)

# Check required number of arguments (at least 7)
if (length(args) < 7) {
  stop("Usage: Rscript script.R <gtf_ref_loc> <gtf_novel_file> <tracking_file> <min_occurrence> <min_tpm> <outfile> <scripts_dir> [gtf_refseq_basename]")
}

# Assign required arguments
gtf_ref_loc     <- args[1]
gtf_novel_file  <- args[2]
tracking_file   <- args[3]
min_occurrence  <- as.numeric(args[4])
min_tpm         <- as.numeric(args[5])
outfile         <- args[6]
scripts_dir     <- args[7]

# Assign optional argument if present
gtf_refseq_basename <- if (length(args) >= 8) args[8] else NULL

# Source additional required functions
#functions_file <- paste0(scripts_dir, "/filter_annotate_functions.R")
functions_file <- "filter_annotate_functions.R"
source(functions_file)


filterGTF <- function( novel_gtf_path, gtf_ref_path, tracking_file, min_occurrence, min_tpm, output_prefix, gtf_refseq_basename) {

  # IMPORT AND PREPARE GTFs------------
  
  # Load reference GTF file for comparison and annotation purposes
  # v2.0: Added exon_number and transcript_biotype, removed transcript_name
  # v2.1: Import GTF allowing either gene_biotype or gene_type to be present
  
  gtf_reference <- rtracklayer::import.gff(gtf_ref_path, colnames = c(
    "type", "source", "gene_id", "gene_name", "gene_biotype", "gene_type", "transcript_id", "transcript_biotype", "exon_number"
  )) 
  
  # Load assembled (novel) GTF file for filtering and fixing
  novel_gtf <- rtracklayer::import(novel_gtf_path)
  
  # Force chromosome style to Ensembl
  GenomeInfoDb::seqlevelsStyle(gtf_reference) <- "Ensembl"
  GenomeInfoDb::seqlevelsStyle(novel_gtf) <- "Ensembl"
  
  # Remove scaffolds and unstranded
  novel_gtf <- novel_gtf[seqnames(novel_gtf) %in% c(1:22, "X", "Y"), ]
  novel_gtf <- novel_gtf[strand(novel_gtf) != "*", ]
  
  # Convert GTF to a data frame for further processing
  gtf_ref_df <- as.data.frame(gtf_reference)
  novel_gtf_df <- data.frame(novel_gtf)
  
  #v2.1: FIX gene_biotype vs gene_type
  
  # If gene_biotype is missing or all NA, but gene_type exists, use gene_type
  if (!"gene_biotype" %in% colnames(gtf_ref_df) ||
      all(is.na(gtf_ref_df$gene_biotype))) {
    if ("gene_type" %in% colnames(gtf_ref_df)) {
      gtf_ref_df$gene_biotype <- gtf_ref_df$gene_type
    }
  }
  
  # Same correction for GRanges metadata
  if (!"gene_biotype" %in% names(mcols(gtf_reference)) ||
      all(is.na(mcols(gtf_reference)$gene_biotype))) {
    if ("gene_type" %in% names(mcols(gtf_reference))) {
      mcols(gtf_reference)$gene_biotype <- mcols(gtf_reference)$gene_type
    }
  }
  
  # Remove gene_type entirely to avoid leaking into output
  gtf_ref_df$gene_type <- NULL
  if ("gene_type" %in% names(mcols(gtf_reference))) {
    mcols(gtf_reference)$gene_type <- NULL
  }
  
  # SELECT NOVEL TRANSCRIPTS ---------------------------
  
  # Discard unwanted transcript classes (reference overlaps) for efficient processing
  unwanted_transctripts <- c("=", "c", "j", "m", "n", "e", "r", "s")
  transcripts_discard <- novel_gtf_df[which(novel_gtf_df$type == "transcript" &
                                              novel_gtf_df$class_code %in% unwanted_transctripts), ]
  novel_gtf_df <- novel_gtf_df[!novel_gtf_df$transcript_id %in% transcripts_discard$transcript_id, ]
  
  #v2.1: Check for empty df after filtering steps
  if (exit_if_empty(novel_gtf_df, gtf_ref_df, output_prefix, script_version, script_date)) {
    return()
  }
  
  # LOGGING: Starting transcript number
  log_start <- length(unique(novel_gtf_df$transcript_id))
  
  # Load tracking file and merge information with novel GTF data frame
  tracking <- read_tracking_file(tracking_file)
  novel_gtf_df <- merge_tracking_info(novel_gtf_df, tracking)
  novel_gtf_df <- fill_metadata_from_transcripts(novel_gtf_df, fields_to_fill = c(
    "gene_name", "oId", "cmp_ref", "class_code", "cmp_ref_gene", "ref_gene_id", "num_samples"
  ))
  novel_gtf_df <- update_gene_id_and_name(novel_gtf_df)
  
  # TRANSCRIPTOME FILTERING ---------------------------
  
  # FILTER: Remove mono-exonic transcripts
  mono_exonic_novel <- count_mono_exonics(gtf = novel_gtf_df)
  novel_gtf_df <- novel_gtf_df[!(novel_gtf_df$transcript_id %in% mono_exonic_novel$transcript_id), ]  
  if (exit_if_empty(novel_gtf_df, gtf_ref_df, output_prefix, script_version, script_date)) {
    return()
  }
  
  # LOGGING: Number of non mono-exonic transcripts
  log_monoexonic <- length(unique(mono_exonic_novel$transcript_id))
  


  # FILTER: Keep transcripts based on minimum occurrence and expression threshold
  # If min_tpm=0 the number of samples from tracking file will be taken
  if(min_tpm > 0) {
    occurrence_mask <- filter_tpm_occurrence(novel_gtf_df, min_occurrence = min_occurrence , min_tpm = min_tpm)
  } else {
    occurrence_mask <- as.numeric(novel_gtf_df$num_samples) > min_occurrence
  }
  
  transcripts_keep <- novel_gtf_df[occurrence_mask, ]$transcript_id
  
  novel_gtf_df <- novel_gtf_df[occurrence_mask, ]
  
  if (exit_if_empty(novel_gtf_df, gtf_ref_df, output_prefix, script_version, script_date)) {
    return()
  }
  
  # LOGGING: Store transcripts failing occurrence threshold
  log_occurrence <- log_start - log_monoexonic - length(unique(transcripts_keep))
  
  # Add biotype to custom annotation based on reference ID
  novel_gtf_df$gene_biotype <- gtf_ref_df$gene_biotype[match(novel_gtf_df$gene_id, gtf_ref_df$gene_id)]
  novel_gtf_df[which(is.na(novel_gtf_df$gene_biotype)), "gene_biotype"] <- "stringtie"
  
  # FILTER: Transcripts overlapping reference transcripts, multiple genes and same strand I class
  unexpected_overlap_txs <- suppressWarnings(find_unexpected_overlaps(gtf = novel_gtf_df, gtf_reference = gtf_reference))
  multigene_txs <- suppressWarnings(find_multigene_overlaps(gtf = novel_gtf_df, gtf_reference = gtf_reference))
  
  ## v2.0: Performance improvement by general findOverlaps() instead of iterating over tx IDs
  same_strand_i_txs <- suppressWarnings(filter_i_class(gtf_df = novel_gtf_df, reference_granges = gtf_reference))
  
  overlapping_txs <- c(unexpected_overlap_txs, multigene_txs, same_strand_i_txs)
  
  novel_gtf_df <- subset(novel_gtf_df,
                         !(novel_gtf_df$transcript_id %in% overlapping_txs))
  
  #v2.1: Check for empty df after filtering steps
  if (exit_if_empty(novel_gtf_df, gtf_ref_df, output_prefix, script_version, script_date)) {
    return()
  }
  
  # LOGGING: Get final number of transcripts to be added
  log_unexpected_overlaps <-  length(unique(unexpected_overlap_txs))
  log_multigene <- length(unique(multigene_txs))
  log_same_strand <- length(unique(same_strand_i_txs))
  log_final <- length(unique(novel_gtf_df$transcript_id))
  
  # WARNINGS: Check for genes with transcripts on multiple strands
  novel_gtf_multiple_strands <- novel_gtf_df %>%
    group_by(gene_id) %>%
    summarize(nstrands = length(unique(strand)))
  
  # Count the number of genes with transcripts on multiple strands
  log_multiple_strands <- sum(novel_gtf_multiple_strands$nstrands > 1)
  
  # WARNINGS: Flagging RefSeq transcripts (only if provided)
  if (!is.null(gtf_refseq_basename) && gtf_refseq_basename != "") {
    novel_gtf_df <- annotate_overlap(gtf = novel_gtf_df, gtf_refseq_basename = gtf_refseq_basename, x_name = "xr")
    novel_gtf_df <- annotate_overlap(gtf = novel_gtf_df, gtf_refseq_basename = gtf_refseq_basename, x_name = "nr")
    novel_gtf_df <- data.frame(novel_gtf_df)
    
    # Optionally report overlap counts
    log_xr_overlap <- nrow(subset(novel_gtf_df, type == "transcript" & xr_overlap == "xr_hit"))
    log_nr_overlap <- nrow(subset(novel_gtf_df, type == "transcript" & nr_overlap == "nr_hit"))
  } else {
    log_xr_overlap <- 0
    log_nr_overlap <- 0
    message("Skipping RefSeq annotation: no gtf_refseq_basename provided.")
  }
  
  # PREPARE OUTPUTS -------------
  ## TODO can this tsv have the same headers as the sample IDs?
  expression_columns <- colnames(novel_gtf_df)[grepl("TPM_q", colnames(novel_gtf_df))]
  output_info <- unique(novel_gtf_df[, c("transcript_id", "gene_id", expression_columns)])
  
  # Define base columns
  base_columns <- c(
    "seqnames", "source", "type", "start", "end", "score", "strand", "phase",
    "gene_id", "transcript_id", "exon_number", "gene_name", "gene_biotype",
    "oId", "cmp_ref", "class_code", "tss_id", "num_samples",
    "cmp_ref_gene", "contained_in", "xloc"
  )
  
  # Add XR/NR columns only if RefSeq path is available
  if(!is.null(gtf_refseq_basename) && gtf_refseq_basename != "") {
    base_columns <- c(base_columns, "xr_overlap", "nr_overlap")
  }
  
  # Subset columns and merge novel with reference
  merged_gtf <- as.data.frame(novel_gtf_df)[, base_columns] %>% 
    bind_rows(gtf_ref_df)
  
  # Sort gtf by feture ranking and coordinate
  # v2.1: Sort both positive and negative strand from 5' to 3'
  merged_gtf_sorted <- sort_gtf(merged_gtf)
  
  # Fill in missing fields
  # v2.1: Add NA string for missing gene/transcript biotypes for compatibility
  merged_gtf_sorted$gene_name <- ifelse(is.na(merged_gtf_sorted$gene_name), merged_gtf_sorted$gene_id,merged_gtf_sorted$gene_name )
  merged_gtf_sorted <- merged_gtf_sorted %>% 
    # Replace missing/empty biotype fields
    mutate(
      gene_name = ifelse(is.na(gene_name),
                         gene_id,
                         gene_name), 
      gene_biotype = if_else(
        is.na(gene_biotype) | gene_biotype == "",
        "NA",
        gene_biotype
      ),
      transcript_biotype = if_else(
        type != "gene" & (is.na(transcript_biotype) | transcript_biotype == ""),
        "NA",
        transcript_biotype
      )
    )
  
  # Turn df into GRanges
  output_gtf <- GenomicRanges::makeGRangesFromDataFrame(merged_gtf_sorted, keep.extra.columns = T)
  
    
  # DEFINE OUTPUT FILES ---------------------------
  output_gtf_path <- paste(output_prefix, "gtf", sep = ".")
  output_info_path <- paste(output_prefix, "tsv", sep = ".")
  output_log_path <- paste(output_prefix, "log", sep = ".")
  
  # WRITE LOG FILE ---------------------------
  
  ## Initialise log file
  cat("#GTF FILTERING", script_version, "(last updated", script_date, ")", "\n", file = output_log_path)
  cat("#gtf_ref_path=", gtf_ref_path, "\n", file = output_log_path, append = T)
  cat("#novel_gtf_path=", novel_gtf_path, "\n", file = output_log_path, append = T)
  cat("#tracking_file=", tracking_file, "\n", file = output_log_path, append = T)
  cat("#gtf_refseq_basename=", gtf_refseq_basename, "\n", file = output_log_path, append = T)
  cat("#min_occurrence=", min_occurrence, "\n", file = output_log_path, append = T)
  cat("#min_tpm=", min_tpm, "\n", file = output_log_path, append = T)
  
  ## Add filtering steps
  cat("Novel transcripts - stranded, not in scaffolds:", "\t", log_start, "\n", file = output_log_path, append = T)
  cat("Filtered monoexonic:", "\t", log_monoexonic, "\n", file = output_log_path, append = T)
  cat("Filtered below expression and occurrence requirement:", "\t", log_occurrence, "\n", file = output_log_path, append = T)
  cat("Filtered reference overlap (class u,x,i,y) :", "\t", log_unexpected_overlaps, "\n", file = output_log_path, append = T)
  cat("Filtered reference multi-gene overlap:", "\t", log_multigene, "\n", file = output_log_path, append = T)
  cat("Filtered class I same strand overlap:", "\t", log_same_strand, "\n", file = output_log_path, append = T)
  cat("__________________________________________", "\n", file = output_log_path, append = T)
  cat("TOTAL NOVEL TRANSCRIPTS ADDED TO REFERENCE:", "\t", log_final, "\n", file = output_log_path, append = T)
  cat("__________________________________________","\n", file = output_log_path, append = T)
  ## Add extra warnings
  cat("Info - Genes with transcripts on multiple strands:", "\t", log_multiple_strands, "\n", file = output_log_path, append = T)
  cat("Info - RefSeq xr overlap:", "\t", log_xr_overlap, "\n", file = output_log_path, append = T)
  cat("Info - RefSeq nr overlap:", "\t", log_nr_overlap, "\n", file = output_log_path, append = T)
  
  # WRITE FINAL GTF AND INFO TABLE ---------------------------
  message(paste(Sys.time(), "Exporting custom gtf ... "), sep = "\t")
  export(object = output_gtf, con = output_gtf_path, format = "gtf", version = "2")
  lines <- readLines(output_gtf_path)
  modified_lines <- gsub("; ID.*", "", lines)
  header_line <- paste0("##Custom GTF generated with filter_annotate.R version ", script_version,
                        " (", script_date, ")")
  modified_lines <- c(header_line, modified_lines)
  writeLines(modified_lines, output_gtf_path)
  write.table(output_info, file = output_info_path, sep = "\t", row.names = F, quote = F)
  message(paste(Sys.time(), "Exporting custom gtf ... Done!"), sep = "\t")
}

# EXECUTE FILTERING  ========================================================
filterGTF(
  novel_gtf_path = gtf_novel_file,
  gtf_ref_path = gtf_ref_loc,
  tracking_file = tracking_file,
  min_occurrence = min_occurrence,
  min_tpm = min_tpm,
  output_prefix = outfile,
  gtf_refseq_basename = gtf_refseq_basename
)
