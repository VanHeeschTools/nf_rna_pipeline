# GTF filtering functions
#' Identify Monoexonic Novel Transcripts from StringTie
#'
#' This function identifies transcripts with only a single exon that were assembled 
#' by StringTie and are marked as novel (i.e., have no reference gene ID).
#' It is typically used to flag transcripts that are likely artefactual or 
#' less informative for downstream analysis.
#'
#' @param gtf A \code{GRanges} object representing the GTF annotation. 
#'   It should include both transcript and exon features, and must contain 
#'   the metadata columns \code{transcript_id}, \code{source}, and \code{ref_gene_id}.
#'
#' @return A \code{GRanges} object containing only the monoexonic transcripts 
#'   that meet the following criteria:
#'   \itemize{
#'     \item Only one exon (i.e., fewer than 2 features per transcript including transcript row)
#'     \item Assembled by StringTie (\code{source == "StringTie"})
#'     \item Novel transcripts (\code{ref_gene_id == "-"})
#'   }
#' 
#' @examples
#' \dontrun{
#' monoexonic <- count_mono_exonics(novel_gtf)
#' }
#'
#' @export
count_mono_exonics <- function(gtf) {
  # Count number of entries per transcript (transcript + exons)
  transcript_counts <- table(gtf$transcript_id) - 1
  
  # Identify transcripts with fewer than 2 entries (monoexonic)
  monoexonic_tx_ids <- names(transcript_counts[transcript_counts < 2])
  
  # Filter GRanges for monoexonic transcripts from StringTie with no reference gene ID
  monoexonic_transcripts <- subset(
    gtf,
    transcript_id %in% monoexonic_tx_ids &
      source == "StringTie" &
      ref_gene_id == "-"
  )
  
  return(monoexonic_transcripts)
}

#' Filter Transcripts by Occurrence and TPM Threshold
#'
#' Filters a data frame of transcript tracking information (e.g., from GFFcompare) 
#' to retain only transcripts that are expressed above a TPM threshold in a minimum 
#' number of samples.
#'
#' @param gtf_df_tracking A \code{data.frame} containing transcript-level information, 
#'   including TPM values per sample. TPM columns should be named starting with \code{"TPM_q"}.
#' @param min_occurrence Integer. Minimum number of samples in which a transcript must be expressed 
#'   at or above \code{min_tpm} to be retained. Default is 1.
#' @param min_tpm Numeric. Minimum TPM required in a sample for it to be counted toward the occurrence. Default is 1.
#'
#' @return A logical vector of the same length as the number of rows in \code{gtf_df_tracking}, 
#'   indicating which transcripts meet the filtering criteria.
#'
#' @examples
#' \dontrun{
#' keep_rows <- filter_tpm_occurrence(gtf_df_tracking, min_occurrence = 2, min_tpm = 1)
#' filtered_df <- gtf_df_tracking[keep_rows, ]
#' }
#'
#' @export
filter_tpm_occurrence <- function(gtf_df_tracking,
                                 min_occurrence = 1,
                                 min_tpm = 1) {
  # Identify TPM columns
  tpm_cols <- grep("^TPM_q", colnames(gtf_df_tracking), value = TRUE)
  
  # Extract TPM matrix and replace NA with 0
  tpm_matrix <- as.matrix(gtf_df_tracking[, ..tpm_cols, drop = FALSE])
  tpm_matrix[is.na(tpm_matrix)] <- 0
  tpm_matrix <- apply(tpm_matrix, 2, as.numeric)
  
  # Count number of samples where TPM >= min_tpm
  tpm_count <- rowSums(tpm_matrix >= min_tpm)
  
  # Return logical vector indicating which transcripts meet the occurrence threshold
  to_keep <- tpm_count >= min_occurrence
  return(to_keep)
}


#' Find Novel Transcripts That Unexpectedly Overlap Reference Exons
#'
#' Identifies transcripts labeled as non-overlapping (class codes 'u', 'x', 'i', 'y')
#' that in fact overlap reference gene exons.
#'
#' @param gtf A data.frame or GRanges-like object with transcript annotations.
#' @param gtf_reference A data.frame or GRanges-like object with reference transcript annotations.
#'
#' @return A character vector of transcript IDs (class codes 'u', 'x', 'i', 'y') overlapping reference exons.
#' 
#' @export
find_unexpected_overlaps <- function(gtf, gtf_reference) {
  gtf_no_overlap <- gtf[gtf$class_code %in% c("u", "x", "i", "y"), ]
  no_overlap_tx <- unique(gtf_no_overlap$transcript_id)
  
  novel_exons <- GenomicRanges::makeGRangesFromDataFrame(
    gtf[gtf$type == "exon" & gtf$transcript_id %in% no_overlap_tx, ],
    keep.extra.columns = TRUE
  )
  ref_exons <- GenomicRanges::makeGRangesFromDataFrame(
    gtf_reference[gtf_reference$type == "exon", ],
    keep.extra.columns = TRUE
  )
  
  hits <- GenomicRanges::findOverlaps(novel_exons, ref_exons, type = "any")
  
  if (length(hits) > 0) {
    overlapped_tx <- unique(mcols(novel_exons)$transcript_id[from(hits)])
      warning(length(overlapped_tx), " transcripts labeled as 'no-overlap' actually overlap reference exons.")
      return(overlapped_tx)
  } else {
    return(character(0))
  }
}


#' Identify Transcripts That Overlap Multiple Reference Genes
#'
#' Finds novel transcripts that overlap exons from more than one reference gene.
#'
#' @param gtf A data.frame or GRanges-like object with transcript annotations.
#' @param gtf_reference A data.frame or GRanges-like object with reference transcript annotations.
#' 
#' @return A character vector of transcript IDs overlapping multiple gene IDs.
#' 
#' @export
find_multigene_overlaps <- function(gtf, gtf_reference) {
  novel_exons <- GenomicRanges::makeGRangesFromDataFrame(
    gtf[gtf$type == "exon", ],
    keep.extra.columns = TRUE
  )
  ref_exons <- GenomicRanges::makeGRangesFromDataFrame(
    gtf_reference[gtf_reference$type == "exon", ],
    keep.extra.columns = TRUE
  )
  
  overlaps <- GenomicRanges::findOverlaps(novel_exons, ref_exons, type = "any")
  if (length(overlaps) == 0) return(character(0))
  
  overlap_df <- data.frame(
    transcript_id = mcols(novel_exons)$transcript_id[from(overlaps)],
    gene_id = mcols(ref_exons)$gene_id[to(overlaps)]
  )
  
  overlap_summary <- split(overlap_df, overlap_df$transcript_id)
  overlap_counts <- vapply(overlap_summary, function(x) length(unique(x$gene_id)), integer(1))
  

  multi_count <- sum(overlap_counts > 1)
  if (multi_count > 0) {
    warning(multi_count, " transcripts overlap exons from more than one reference gene.")
    return(names(overlap_counts[overlap_counts > 1]))
  } else {
    return(character(0))
  }
}

#' Filter same-strand I-class transcripts overlapping known reference transcripts
#'
#' This function identifies and returns transcript IDs from a novel transcript GTF
#' that are classified as "i" (fully intronic) and overlap any known reference 
#' transcript on the same strand. Used to filter out likely annotation artifacts.
#'
#' @param gtf_df A data frame containing novel transcripts, including `transcript_id`,
#' `type`, and `class_code` columns.
#' @param reference_granges A GRanges object containing reference transcript annotations.
#'
#' @return A character vector of transcript IDs that overlap known transcripts and 
#' should be removed.
#'
#' @import GenomicRanges
filter_i_class <- function(gtf_df, reference_granges) {
  # Subset for "i" class transcripts
  gtf_i <- subset(gtf_df, type == "transcript" & class_code == "i")
  
  if (nrow(gtf_i) == 0) {
    return(character(0))
  }
  
  # Convert to GRanges
  gtf_i_granges <- GenomicRanges::makeGRangesFromDataFrame(gtf_i, keep.extra.columns = TRUE)
  
  # Find overlaps once for all "i" class transcripts
  hits <- GenomicRanges::findOverlaps(query = gtf_i_granges, subject = reference_granges, type = "any", ignore.strand = FALSE)
  
  if (length(hits) == 0) {
    return(character(0))
  } else {
    # Extract the transcript_ids for overlapping transcripts
    overlapping_tx_ids <- gtf_i_granges$transcript_id[queryHits(hits)]
    
    # Return unique transcript_ids
    return(unique(overlapping_tx_ids))
  }

}


#' Annotate Novel Transcripts with Overlaps in RefSeq GFF
#'
#' This function annotates novel transcripts with any overlap in the given RefSeq GFF file.
#' It is specifically used to identify and annotate XR_### and NR_### transcripts that might not 
#' be applicable for neo-antigen detection based on overlaps with known exons.
#'
#' @param gtf A \code{GRanges}-like object representing novel transcript annotations, 
#'   containing exon coordinates and metadata (e.g., from a GTF file).
#' @param gtf_refseq_basename The basename of the RefSeq GFF file (without extension) 
#'   that contains the exon annotations for XR and NR transcripts.
#' @param x_name A character vector containing two elements: 
#'   - The name of the new column to be added to the \code{gtf} object.
#'   - The name of the RefSeq GFF file that contains the transcript annotations to check for overlap.
#'
#' @return A \code{GRanges} object with the novel transcripts annotated with whether or not 
#'   they overlap any exons in the provided RefSeq GFF file.
#'   The new column in the \code{GRanges} object will be named according to \code{x_name},
#'   and it will contain "hit" for overlapping transcripts or "none" for those without any overlap.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Imports the RefSeq GFF file.
#'   \item Filters the RefSeq GFF file based on valid chromosome names (1:22, X, Y).
#'   \item Subsets the RefSeq GFF file to only include exon-type entries.
#'   \item Checks for overlaps between the novel transcripts and the RefSeq exons.
#'   \item Annotates the novel transcripts with a "hit" if they overlap any RefSeq exons, or "none" if they do not.
#' }
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom rtracklayer import
#' @importFrom stats unique
#' @export
annotate_overlap <- function(gtf, gtf_refseq_basename, x_name) {
  x <- rtracklayer::import(paste(gtf_refseq_basename, x_name, "gff", sep = "."))
  
  seqnames_list <- unique(seqnames(x))
  valid_chromosomes <- c(1:22, "X", "Y")
  
  if (any(valid_chromosomes %in% seqnames_list)) {
    x <- x[seqnames(x) %in% valid_chromosomes]
  } else if (any(grepl("^chr[0-9XY]+$", seqnames_list))) {
    seqlevels(x) <- gsub("^chr", "", seqlevels(x))
    x <- x[seqnames(x) %in% valid_chromosomes]
  } else if (any(grepl("^NC_", seqnames_list))) {
    refseq_to_ensembl <- c(
      "NC_000001.11" = "1", "NC_000002.12" = "2", "NC_000003.12" = "3",
      "NC_000004.12" = "4", "NC_000005.10" = "5", "NC_000006.12" = "6",
      "NC_000007.14" = "7", "NC_000008.11" = "8", "NC_000009.12" = "9",
      "NC_000010.11" = "10", "NC_000011.10" = "11", "NC_000012.12" = "12",
      "NC_000013.11" = "13", "NC_000014.9"  = "14", "NC_000015.10" = "15",
      "NC_000016.10" = "16", "NC_000017.11" = "17", "NC_000018.10" = "18",
      "NC_000019.10" = "19", "NC_000020.11" = "20", "NC_000021.9"  = "21",
      "NC_000022.11" = "22", "NC_000023.11" = "X", "NC_000024.10" = "Y"
    )
    x <- x[seqnames(x) %in% names(refseq_to_ensembl)]
    x <- keepSeqlevels(x, value = names(refseq_to_ensembl))
    seqlevels(x) <- refseq_to_ensembl[seqlevels(x)]
  } else {
    stop("Unknown seqname format. Please check the input.")
  }
  
  x <- subset(x, x$type == "exon")
  
  gtf_novel_gr <- GenomicRanges::makeGRangesFromDataFrame(gtf, keep.extra.columns = TRUE)
  gtf_novel_exons <- subset(gtf_novel_gr, gtf_novel_gr$type == "exon")
  
  x_overlap <- GenomicRanges::findOverlaps(query = x, subject = gtf_novel_exons, type = "any")
  
  novel_x_hits <- as.data.frame(unique(gtf_novel_exons[subjectHits(x_overlap)]))
  novel_x_hits <- unique(novel_x_hits[, "transcript_id"])
  
  elementMetadata(gtf_novel_gr)[[paste0(x_name, "_overlap")]] <- ifelse(
    gtf_novel_gr$transcript_id %in% novel_x_hits, paste0(x_name, "_hit"), "none"
  )
  
  return(gtf_novel_gr)
}

#' Rename Transcripts Without Reference Gene
#'
#' This function renames transcripts that do not have a reference gene by generating a new gene name
#' based on the chromosome, start, and end positions of the transcript. Linked transcripts (sharing the same gene) 
#' will have the same gene name.
#'
#' @param gtf_novel_df A data frame or data table containing novel transcript annotations, including 
#'   chromosome, start, and end positions, as well as transcript IDs.
#'
#' @return A data frame with updated gene names for transcripts that do not have a reference gene ID.
#'   The new gene name is a combination of the chromosome, start, and end position.
#'
#' @details This function performs the following steps:
#' \itemize{
#'   \item Splits the novel transcripts by gene.
#'   \item For each gene, if there is no reference gene ID, a new gene name is generated based on 
#'   the chromosome and positions of the transcript.
#'   \item Returns the updated novel transcript annotations with the renamed genes.
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom stats na.omit
#' @export
rename_stringtie_transcripts <- function(gtf_novel_df) {
  gtf_by_gene <- split(gtf_novel_df, gtf_novel_df$gene_id)
  
  gtf_renamed <- suppressWarnings(bind_rows(lapply(gtf_by_gene, function(x) {
    if (is.na(x$ref_gene_id)[1]) {
      chr <- unique(x$seqnames)
      start <- min(x$start)
      end <- max(x$end)
      x$gene_name <- paste0(chr, ":", start, "-", end)
    }
    
    return(x)
  })))
}

#' Read and parse a tracking file with transcript annotations and TPMs
#'
#' This function reads a tracking file output (e.g., from StringTie), extracts transcript and gene
#' annotations, and parses transcript-level TPM values from each sample.
#'
#' @param tracking_file A character string specifying the path to the tracking file (typically a `.txt` file).
#'
#' @return A data frame containing:
#' \describe{
#'   \item{TCONS}{Transcript ID (e.g., TCONS_xxxx)}
#'   \item{xloc}{Locus identifier (e.g., XLOC_xxxx)}
#'   \item{ref_gene_id}{Reference gene ID parsed from column 3}
#'   \item{ref_transcript_id}{Reference transcript ID parsed from column 3}
#'   \item{TPM_qX}{TPM values extracted from sample columns, one column per sample (e.g., `TPM_q1`, `TPM_q2`, ...)}
#' }
#'
#' @details The function expects the third column of the file to contain `ref_gene_id|ref_transcript_id` format.
#' Remaining columns are parsed to extract TPM values from the third-to-last field in a `|`-separated string.
#'
#' Missing or malformed values (`"-"` or too few fields) are returned as `NA`.
#'
#' @importFrom tidyr separate
#' @importFrom dplyr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df <- read_tracking_file("stringtie/merge_tracking.txt")
#' head(df)
#' }
read_tracking_file <- function(tracking_file) {
  # Read the tracking file
  df <- read.table(tracking_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract the first 4 columns
  df_cleaned <- df[, 1:3] %>%
    tidyr::separate(V3,
                    into = c("ref_gene_id_annotated", "ref_transcript_id_annotated"),
                    sep = "\\|"
    )
  colnames(df_cleaned) <- c("TCONS", "xloc", "ref_gene_id", "ref_transcript_id")
  
  # Extract TPMs from the remaining columns
  df_tpm <- df[, 5:ncol(df)] %>%
    apply(1, function(row) {
      sapply(row, function(x) {
        if (x == "-" || is.na(x)) {
          return(NA)
        }
        parts <- unlist(strsplit(x, "\\|"))
        if (length(parts) >= 3) {
          return(parts[length(parts) - 2])
        } else {
          return(NA)
        }
      })
    }) %>%
    t() %>%
    as.data.frame()
  
  # Add sample names (q1, q2, ...) as column names for the TPMs
  colnames(df_tpm) <- paste0("TPM_q", 1:ncol(df_tpm))
  
  # Combine the first 4 columns and the extracted TPMs into a final data frame
  final_df <- cbind(df_cleaned, df_tpm)
  
  return(final_df)
}

#' Merge tracking info into a GTF data frame
#'
#' @param gtf_df A data frame version of the novel GTF.
#' @param tracking_df A data frame parsed from a tracking file.
#'
#' @return A data frame with tracking columns merged in by transcript_id.
merge_tracking_info <- function(gtf_df, tracking_df) {
  dplyr::left_join(gtf_df, tracking_df, by = c("transcript_id" = "TCONS"))
}

#' Fill missing metadata in exon rows with values from transcript rows
#'
#' @param gtf_df A data.table with transcript and exon rows.
#' @param fields_to_fill A character vector of columns to carry over from transcript rows.
#'
#' @return The input data.table with NA values filled.
fill_metadata_from_transcripts <- function(gtf_df, fields_to_fill) {
  data.table::setDT(gtf_df)
  gtf_df[, (fields_to_fill) := {
    w <- which(type == "transcript")
    if (length(w) == 0) return(.SD)
    lapply(.SD, function(x) {
      x[is.na(x)] <- x[w]
      x
    })
  }, by = transcript_id, .SDcols = fields_to_fill]
  return(gtf_df)
}

#' Update gene_id and gene_name using reference gene info
#'
#' @param gtf_df A data.frame with columns `ref_gene_id`, `cmp_ref_gene`, `gene_id`, and `gene_name`.
#'
#' @return The input data.frame with updated `gene_id` and `gene_name`.
update_gene_id_and_name <- function(gtf_df) {
  gtf_df$gene_id <- ifelse(
    gtf_df$ref_gene_id != "-" & !is.na(gtf_df$ref_gene_id) & is.na(gtf_df$cmp_ref_gene) ,
    gtf_df$ref_gene_id,
    gtf_df$gene_id
  )
  gtf_df$gene_name <- ifelse(
    is.na(gtf_df$gene_name),
    gtf_df$gene_id,
    gtf_df$gene_name
  )
  return(gtf_df)
}

#' Sort GTF Entries by Gene, Transcript, and Feature Hierarchy
#'
#' Sorts a GTF data.frame or data.table by genomic location and feature hierarchy:
#' 1. Genes globally by chromosome, start, and reverse end.
#' 2. Transcripts within each gene by start and reverse end.
#' 3. Features (exon/CDS/UTR/others) within each transcript, with known types first,
#'    and strand-aware ordering by start (or reversed start for negative strand).
#'
#' @param gtf_df A data.frame or data.table containing GTF fields: seqnames, start, end,
#'   strand, gene_id, transcript_id, and type (feature type).
#' @return A data.table with the same columns as input, sorted according to the hierarchy.
#' @export
sort_gtf <- function(gtf_df) {
  #Transform to DT for improved speed
  dt <- data.table::as.data.table(gtf_df)
  
  # Step 1: Rank genes by coordinate
  # Get existing gene entries
  gene_rows <- dt[type == "gene", .(seqnames, start, end, strand, gene_id)]
  
  # Identify missing gene_ids
  all_gene_ids <- unique(dt$gene_id)
  gene_ids_missing_row <- setdiff(all_gene_ids, gene_rows$gene_id)
  
  # Synthesize gene rows from representative transcripts
  synth_gene_rows <- dt[
    gene_id %in% gene_ids_missing_row & type == "transcript",
    .SD[which.min(start)],  # representative transcript per gene
    by = gene_id
  ][, .(seqnames, start, end, strand, gene_id)]
  
  # Combine for ranking
  gene_dt <- rbind(gene_rows, synth_gene_rows)
  data.table::setorder(gene_dt, seqnames, start, -end)
  gene_dt[, gene_rank := .I]
  
  # Merge gene_rank
  dt <- gene_dt[, .(gene_id, gene_rank)][dt, on = "gene_id"]
  
  # Step 2: Rank transcripts within genes
  trans_dt <- dt[type == "transcript", .(gene_id, transcript_id,seqnames, start, end)]
  trans_dt <- unique(trans_dt)
  data.table::setorder(trans_dt, seqnames,  start, -end, gene_id, transcript_id)
  trans_dt[, transcript_rank := seq_len(.N), by = gene_id]
  
  # Merge transcript_rank
  dt <- trans_dt[, .(transcript_id, transcript_rank)][dt, on = "transcript_id"]
  
  # Step 3: Feature type priority
  priority_map <- c(
    gene = 1,
    transcript = 2,
    exon = 3,
    CDS = 4,
    UTR = 4, 
    start_codon = 4
  )
  
  dt[, feature_rank := data.table::fifelse(
    type %in% names(priority_map),
    priority_map[as.character(type)],
    4 + as.integer(as.factor(type))
  )]
  
  # Step 4: Start-aware feature position
  dt[type %in% c("exon", "CDS", "UTR", "start_codon") & strand == "+",  feature_order := start]
  dt[type %in% c("exon", "CDS", "UTR", "start_codon") & strand == "-",  feature_order := -start]
  
  dt[!type %in% c("exon", "CDS", "UTR", "start_codon"), feature_order := start]
  
  # Step 5: Final sort
  data.table::setorder(dt, gene_rank, transcript_rank, feature_rank, feature_order)
  
  # Step 6: Cleanup temporary columns
  dt[, c("gene_rank", "transcript_rank", "feature_rank", "feature_order") := NULL]
  
  return(dt)
}


#' Exit Early if Novel GTF is Empty and Write Reference-Only Output
#'
#' This helper function checks whether the `novel_gtf_df` has become empty during
#' the filtering pipeline. If it is empty, the function:
#' \itemize{
#'   \item Emits a message informing that only the reference GTF will be exported.
#'   \item Exports the reference GTF to the specified output path.
#'   \item Prepends a header line with script metadata (version and date).
#'   \item Writes an empty info table to the specified output path.
#'   \item Writes a log file describing that the novel entries became empty.
#'   \item Returns `TRUE` to indicate that execution should stop.
#' }
#'
#' If `novel_gtf_df` is not empty, the function returns `FALSE`, allowing the
#' pipeline to continue normally.
#'
#' @param novel_gtf_df A data frame containing novel GTF entries.
#' @param gtf_ref_df A data frame containing the reference GTF to be exported
#'        if no novel entries remain.
#' @param output_gtf_path File path where the resulting GTF should be written.
#' @param output_log_path File path for writing the log file describing the exit event.
#' @param script_version A character string indicating the script version.
#' @param script_date A character string indicating the script date.
#'
#' @return `TRUE` if the function exported the reference GTF and the caller should
#'         exit further processing; `FALSE` otherwise.
#'
#' @import GenomicRanges
#' @importFrom rtracklayer export
#' @export
exit_if_empty <- function(novel_gtf_df, gtf_ref_df, output_prefix, script_version, script_date) {
  # DEFINE OUTPUT FILES ---------------------------
  output_gtf_path <- paste(output_prefix, "gtf", sep = ".")
  output_info_path <- paste(output_prefix, "tsv", sep = ".")
  output_log_path <- paste(output_prefix, "log", sep = ".")
  
  if (nrow(novel_gtf_df) == 0) {
    
    message("novel_gtf_df is empty - exporting reference-only output.")
    
    # ---- Write reference GTF with header ----
    output_gtf <- GenomicRanges::makeGRangesFromDataFrame(
      gtf_ref_df, keep.extra.columns = TRUE
    )
    export(output_gtf, con = output_gtf_path, format = "gtf", version = "2")
    
    lines <- readLines(output_gtf_path)
    writeLines(
      c(
        paste0("##GTF generated with filter_annotate.R version ",
               script_version, " (", script_date, ")"),
        lines
      ),
      output_gtf_path
    )
    
    # ---- Write empty info table ----
    write.table(
      data.frame(),
      file = output_info_path,
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    
    # ---- Write log ----
    cat(
      "#GTF FILTERING", script_version, "(updated ", script_date, ")\n",
      "Novel GTF became empty during filtering.\n",
      "Output contains only reference GTF.\n",
      file = output_log_path,
      sep = ""
    )
    
    return(TRUE)
  }
  
  return(FALSE)
}
