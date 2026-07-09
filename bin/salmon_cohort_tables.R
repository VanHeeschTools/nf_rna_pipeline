#!/usr/bin/env Rscript

# Author: Marina Reixachs Sole [M.Reixachssole@prinsesmaximacentrum.nl]
# Usage: script.R <quant_paths> <gtf_file> <prefix>

suppressPackageStartupMessages({
  library(rtracklayer)
  library(tximport)
  library(readr)
  library(dplyr)
  library(stringr)
})

# READ ARGS-----------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: script.R <quant_paths> <gtf_file> <prefix>")
}

quant_paths <- args[1]
gtf_file <- args[2]
prefix <- args[3]

print(quant_paths)

# LOGGING FUNCTION -----------------
log_file <- file.path(paste0(prefix, "_log.txt"))
log_message <- function(message, level = "WARNING") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", level, "] ", timestamp, " - ", message)
  write(log_entry, file = log_file, append = TRUE)
  message(log_entry) # Also print to console for visibility
}

# FUNCTIONS -----------------

# Function to generate tx2gene
get_tx2gene <- function(gtf_data) {
  # Select transcript entries only
  transcript_entries <- as.data.frame(gtf_data) %>%
    filter(type == "transcript")
  
  # Stop if no transcript entries are found
  if (is.null(transcript_entries) || nrow(transcript_entries) == 0) {
    stop("No valid tx2gene data. Ensure the GTF file contains transcript entries.")
  }
  
  # Evaluate presence of gene_id and gene_name
  has_gene_id <- "gene_id" %in% colnames(transcript_entries) && !all(is.na(transcript_entries$gene_id))
  has_gene_name <- "gene_name" %in% colnames(transcript_entries) && !all(is.na(transcript_entries$gene_name))
  
  # Error handling for missing fields
  if (!has_gene_id && !has_gene_name) {
    log_message("No 'gene_id' or 'gene_name' found in the GTF file. Gene-level expression will not be computed.")
    return(NULL)
  }
  
  # Select and handle gene_id and gene_name
  if (all(is.na(transcript_entries$gene_name))) {
    log_message("All 'gene_name' values are NA. Gene name table will not be generated.")
    tx2gene <- tx2gene %>% 
      dplyr::select(transcript_id, gene_id) %>%
      dplyr::mutate(transcript_id = gsub("\\..*", "", transcript_id),
                    gene_id = gsub("\\..*", "", gene_id))
  } else if (all(is.na(transcript_entries$gene_id))) {
    log_message("All 'gene_id' values are NA. Gene ID table will not be generated.")
    tx2gene <- transcript_entries %>%
      dplyr::select(transcript_id, gene_name) %>%
      dplyr::mutate(transcript_id = gsub("\\..*", "", transcript_id))
  } else {
    tx2gene <- transcript_entries %>%
      dplyr::select(transcript_id, gene_id, gene_name) %>%
      dplyr::mutate(transcript_id = gsub("\\..*", "", transcript_id),
                    gene_id = gsub("\\..*", "", gene_id)) %>%
      mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%  # Replace NA gene_name with gene_id
      mutate(gene_id = ifelse(is.na(gene_id), gene_name, gene_id)) 
  }
  
  # Remove rows with missing transcript IDs
  tx2gene <- tx2gene[!is.na(tx2gene$transcript_id), ]
  
  # Log missing values
  missing_gene_ids <- sum(is.na(tx2gene$gene_id))
  if (missing_gene_ids > 0) {
    log_message(paste(missing_gene_ids, "transcripts have missing 'gene_id' values."))
  }
  
  if ("gene_name" %in% colnames(tx2gene)) {
    missing_gene_names <- sum(is.na(tx2gene$gene_name))
    if (missing_gene_names > 0) {
      log_message(paste(missing_gene_names, "transcripts have missing 'gene_name' values."))
    }
  }
  
  return(tx2gene)
}
# Function to write properly formatted tables
write_tsv <- function(data, file) {
  if (!is.null(data)) {
    write.table(data,
                file = file,
                sep = "\t",
                row.names = TRUE,
                col.names = TRUE,
                quote = FALSE
    )
    log_message(paste("Successfully wrote file:", file), level = "INFO")
  } else {
    log_message(paste("Unable to write output file:", file))
  }
}

# Function to import salmon files - aware of gencode transcript IDs (ENST..|ENSG..)
import_salmon <- function(quant_files) {
  # Load first row to check format in salmon
  first_record <- read_tsv(quant_files[1], n_max = 1, show_col_types = FALSE)
  first_id <- first_record$Name
  
  if (grepl("\\|", first_id)) {
    txi <- tximport(
      files = quant_files,
      type = "salmon",
      txOut = TRUE,
      ignoreAfterBar = TRUE,
      dropInfReps = TRUE
    )
   
  } else {
    txi <- tximport(
      files = quant_files,
      type = "salmon",
      txOut = TRUE,
      dropInfReps = TRUE
    )
  }
  
  #Manually remove tx version (ignoreAfterBar and ignoreTxVersion are incompatible)
  rownames(txi$abundance) <- gsub("\\..*", "", rownames(txi$abundance)) 
  rownames(txi$length) <- gsub("\\..*", "", rownames(txi$length)) 
  rownames(txi$counts) <- gsub("\\..*", "", rownames(txi$counts)) 
  
  if (is.null(txi$counts) || is.null(txi$abundance)) {
    stop("Failed to load quantification data. Check your Salmon outputs and paths.")
  }
  
  return(txi)
}

# CODE -----------------

# Create or clear log file
write("", file = log_file)

# Import gtf_file
gtf_data <- rtracklayer::import(gtf_file)

# Get tx2gene object
tx2gene <- get_tx2gene(gtf_data)

# Save tx2gene table to file
tx2gene_file <- file.path(paste0(prefix, "_tx2gene.tsv"))
write_tsv(tx2gene, tx2gene_file)

quant_files <- readLines(quant_paths)

# Get the list of quant.sf files
if (length(quant_files) == 0) {
  log_message("No quant.sf files found. Ensure the input file contains valid paths.", level = "ERROR")
  stop("No quant.sf files found. Check your input.")
}

# Create a named vector for files
names(quant_files) <- basename(dirname(quant_files))

# Import the quantification files using tximport
txi <- import_salmon(quant_files)

# Summarise gene counts by gene_id
if ("gene_id" %in% colnames(tx2gene)) {
  if (all(tx2gene$transcript_id == tx2gene$gene_id, na.rm = TRUE)) {
    log_message("All 'transcript_id' values are identical to 'gene_id'. Gene ID tables will not be generated.")
    txi.id <- NULL
  } else {
    txi.id <- summarizeToGene(txi, tx2gene[, c("transcript_id", "gene_id")])
  }
} else {
  txi.id <- NULL
}

# Summarise gene counts by gene_name
if ("gene_name" %in% colnames(tx2gene)) {
  if (all(tx2gene$transcript_id == tx2gene$gene_name, na.rm = TRUE)) {
    log_message("All 'transcript_id' values are identical to 'gene_name'. Gene name tables will not be generated.")
    txi.name <- NULL
  } else {
    txi.name <- summarizeToGene(txi, tx2gene[, c("transcript_id", "gene_name")])
  }
} else {
  txi.name <- NULL
}

# Define output file names
counts_file <- file.path(paste0(prefix, "_transcript_counts.tsv"))
tpms_file <- file.path(paste0(prefix, "_transcript_tpms.tsv"))
gene_ID_counts_file <- file.path(paste0(prefix, "_gene_ID_counts.tsv"))
gene_ID_tpms_file <- file.path(paste0(prefix, "_gene_ID_tpms.tsv"))
gene_name_counts_file <- file.path(paste0(prefix, "_gene_name_counts.tsv"))
gene_name_tpms_file <- file.path(paste0(prefix, "_gene_name_tpms.tsv"))

# Write outputs
write_tsv(txi$counts, counts_file)
write_tsv(txi$abundance, tpms_file)
write_tsv(txi.id$counts, gene_ID_counts_file)
write_tsv(txi.id$abundance, gene_ID_tpms_file)
write_tsv(txi.name$counts, gene_name_counts_file)
write_tsv(txi.name$abundance, gene_name_tpms_file)


# Multiqc tables

# Per-sample summary (mean and optional median, with and without zeros)
sample_summary <- function(mat, prefix, median = FALSE) {
  df <- data.frame(sample = colnames(mat),
                  #mean_all = colMeans(mat, na.rm = TRUE),
                  mean_bigger_than_1 = apply(mat, 2, function(x) mean(x[x>1], na.rm=TRUE)))
  if (median) {
    df$median_bigger_than_1 <- apply(mat, 2, function(x) median(x[x>1], na.rm=TRUE))
  }
  setNames(df, c("sample", paste0(prefix, "_", names(df)[-1])))
}

# Parse Salmon meta_info.json for multiple samples
salmon_qc <- function(quant_files) {
  cols <- c("num_processed", "num_mapped", "percent_mapped")
  
  do.call(rbind, lapply(names(quant_files), function(smp) {
    f <- file.path(dirname(quant_files[smp]), "aux_info/meta_info.json")
    qc <- if (file.exists(f)) tryCatch(rjson::fromJSON(file = f), error = function(e) list()) else list()
    data.frame(sample = smp, t(vapply(cols, function(x) qc[[x]] %||% NA, numeric(1))),
              stringsAsFactors = FALSE)
  }))
}

multiqc_table <- Reduce(function(x, y) dplyr::full_join(x, y, by="sample"),
                        list(
                          #sample_summary(txi$counts, "transcript_counts"),
                          sample_summary(txi$abundance, "transcript_tpms", median=TRUE),
                          #sample_summary(txi.id$counts, "gene_ID_counts"),
                          sample_summary(txi.id$abundance, "gene_ID_tpms", median=TRUE)
                        ))
multiqc_table$expressed_transcripts <- colSums(txi$abundance > 1)
multiqc_table$expressed_genes <- colSums(txi.id$abundance > 1)
multiqc_table <- dplyr::relocate(
  multiqc_table,
  expressed_transcripts, expressed_genes,
  .before = transcript_tpms_mean_bigger_than_1)

multiqc_table <- dplyr::left_join(salmon_qc(quant_files), multiqc_table, by="sample")

name_map <- c(
  sample = "Sample",
  num_processed = "Reads Processed",
  num_mapped = "Reads Mapped",
  percent_mapped = "Percent Mapped",
  expressed_transcripts = "Transcripts with TPM >1",
  expressed_genes = "Genes with TPM >1",
  #transcript_counts_mean_all = "Transcript Mean Counts",
  #transcript_counts_mean_bigger_than_1 = "Transcript Mean Counts >1",
  #transcript_tpms_mean_all = "Transcript Mean TPM",
  transcript_tpms_mean_bigger_than_1= "Transcript Mean TPM >1",
  transcript_tpms_median_bigger_than_1 = "Transcript Median TPM >1",
  #gene_ID_counts_mean_all = "Gene Mean Counts",
  #gene_ID_counts_mean_bigger_than_1 = "Gene Mean Counts >1",
  #gene_ID_tpms_mean_all = "Gene Mean TPM",
  gene_ID_tpms_mean_bigger_than_1 = "Gene Mean TPM >1",
  gene_ID_tpms_median_bigger_than_1 = "Gene Median TPM >1"
)

multiqc_table <- multiqc_table %>%
  rename_with(~ name_map[.x]) 

write.table(multiqc_table,
            file = paste0(prefix, "_multiqc_summary_mqc.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)



