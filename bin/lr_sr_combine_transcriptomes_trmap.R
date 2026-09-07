#!/usr/bin/env Rscript

# Global arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign required arguments
sr_gtf_file    <- args[1]
lr_gtf_file    <- args[2]
trmap_out_file <- args[3]
trmap_tab_file <- args[4]
qc_output      <- args[5]
min_occurrence <- args[6]

# Generate header early (used consistently throughout)
gtf_header <-  c(
  paste("##UNIFIED LR/SR GTF", 1, "(last updated 30-04-2026)" ),
  paste0("##SR GTF=", sr_gtf_file),
  paste0("##LR GTF=", lr_gtf_file),
  paste0("##trmap file=", trmap_out_file),
  paste0("##min_occurrence=", min_occurrence),
  paste0("##min_tpm_lr=", 0.1),
  paste0("##min_tpm_sr=", 1)
)

# LIBRARIES ------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

functions_file <- "filter_annotate_functions.R"
source(functions_file)

# ── LOAD & PREPRICESS ──────────────────────────────────────
# Load assembled (novel) GTF file for filtering and fixing
sr_gtf <- rtracklayer::import(sr_gtf_file)
lr_gtf <- rtracklayer::import(lr_gtf_file)

sr_gtf_df <- as.data.frame(sr_gtf) %>% rename_with(~ paste(., "sr", sep = "."))
lr_gtf_df <-  as.data.frame(lr_gtf)
mono_exonic_lr <- count_mono_exonics(gtf = lr_gtf_df) %>% pull(transcript_id)
lr_gtf_df <-  lr_gtf_df %>% rename_with(~ paste(., "lr", sep = ".")) %>%  filter(!transcript_id.lr %in% mono_exonic_lr)

read_trmap_out_full <- function(trmap_out_file) {
  
  lines <- readLines(trmap_out_file)
  
  result <- list()
  current_query <- NULL
  
  for (ln in lines) {
    
    # HEADER LINE ----------------------------------
    if (grepl("^>", ln)) {
      
      header <- sub("^>", "", ln)
      fields <- strsplit(header, " ")[[1]]
      
      tx_id   <- fields[1]
      genomic <- fields[2]
      strand  <- fields[3]
      exons   <- fields[4]
      
      chr <- sub(":.*", "", genomic)
      coords <- sub(".*:", "", genomic)
      start <- as.numeric(sub("-.*", "", coords))
      end   <- as.numeric(sub(".*-", "", coords))
      
      current_query <- tibble(
        transcript_id.lr = tx_id,
        seqnames.lr = chr,
        start.lr = start,
        end.lr = end,
        strand.lr = strand,
        exon_chain.lr = exons
      )
      
    } else {
      
      # OVERLAP LINE -------------------------------
      # Format:
      # class  chr  strand  start  end  refID  exon_list
      
      fields <- strsplit(ln, "\t")[[1]]
      
      overlap_row <- tibble(
        cmp_class = fields[1],
        seqnames.sr = fields[2],
        strand.sr = fields[3],
        start.sr = as.numeric(fields[4]),
        end.sr   = as.numeric(fields[5]),
        transcript_id.sr    = as.character(fields[6]),
        exon_chain.sr = fields[7]
      )
      
      # Combine header + overlap
      result[[length(result) + 1]] <- bind_cols(current_query, overlap_row)
    }
  }
  
  bind_rows(result)
}

trmap_df <- read_trmap_out_full(trmap_out_file) 
trmap_tab_df <- read.delim(trmap_tab_file, header = F, stringsAsFactors = F) 
colnames(trmap_tab_df) <- c("transcript_id.lr", "cmp_class", "matched_percent","transcript_id.sr", "matched_junctions")

# Extract reference for later
reference <- sr_gtf_df[sr_gtf_df$source.sr != "StringTie",]
colnames(reference) <- gsub(".sr", "", colnames(reference))

# ── MERGE WITH UNMATCHED ──────────────────────────────────────
all_transcripts <- trmap_tab_df %>%
filter(!transcript_id.lr %in% mono_exonic_lr)   %>% 
  left_join(lr_gtf_df %>%
              filter(source.lr == "StringTie", 
                      type.lr == "transcript") %>% 
              dplyr::select( gene_id.lr, 
                      gene_name.lr, 
                      transcript_id.lr, 
                      num_samples.lr,
                      num_samples_TPM_threshold.lr)) %>%
  left_join(sr_gtf_df %>%
              filter(source.sr == "StringTie", 
                      type.sr == "transcript") %>% 
              dplyr::select( gene_id.sr, 
                      gene_name.sr, 
                      transcript_id.sr, 
                      num_samples.sr,
                      num_samples_TPM_threshold.sr) ) %>%
  bind_rows(lr_gtf_df %>%
              filter(source.lr == "StringTie", 
                      type.lr == "transcript", 
                      !transcript_id.lr %in% trmap_tab_df$transcript_id.lr) %>% 
              dplyr::select( gene_id.lr, 
                      gene_name.lr, 
                      transcript_id.lr, 
                      num_samples.lr,
                      num_samples_TPM_threshold.lr)) %>%
  bind_rows(sr_gtf_df %>%
              filter(source.sr == "StringTie", 
                      type.sr == "transcript", 
                      !transcript_id.sr %in% trmap_tab_df$transcript_id.sr) %>% 
              dplyr::select( gene_id.sr, 
                      gene_name.sr, 
                      transcript_id.sr, 
                      num_samples.sr,
                      num_samples_TPM_threshold.sr))

# ── ASSIGN CLASS PRIORITY ──────────────────────────────────────

trmap_df <- trmap_df %>%
  left_join(trmap_tab_df) %>% 
  filter(!transcript_id.lr %in% mono_exonic_lr) %>%
  left_join(
    lr_gtf_df %>%
      filter(type.lr == "transcript") %>%
      dplyr::select(transcript_id.lr, gene_id.lr, gene_name.lr, num_samples.lr, num_samples_TPM_threshold.lr)) %>%
  left_join(
    sr_gtf_df %>%
      filter(type.sr == "transcript") %>%
      dplyr::select(transcript_id.sr, gene_id.sr, gene_name.sr,  num_samples.sr, num_samples_TPM_threshold.sr),
  ) %>%
  mutate(num_samples_TPM_threshold.lr = as.numeric(num_samples_TPM_threshold.lr),
        num_samples_TPM_threshold.sr = as.numeric(num_samples_TPM_threshold.sr),
        combined_threshold = num_samples_TPM_threshold.lr + num_samples_TPM_threshold.sr,
        recurrence_pass = ifelse(cmp_class == "=",
                                  combined_threshold >= min_occurrence,
                                  (num_samples_TPM_threshold.lr >= min_occurrence) & (num_samples_TPM_threshold.sr >= min_occurrence))
  )

trmap_df_filt <- trmap_df %>%
  filter(recurrence_pass)

# DEPRECATED
# class_priority_rank <- c(
#   "="  = 1, # Complete match, identical intron chain, same TSS and TES
#   "c"  = 2, # Contained within a reference transcript, same splice sites, subset of introns
#   "j"  = 3, # Multi-exon, shares at least one splice junction with a known transcript
#   "k"  = 4, # Containment reverse OR hard to classify cases in the same strand
#   "m"  = 5, # Retained intron, monoexon fragment overlapping an intron
#   "n"  = 5, # Intron retention confirmed by intron chain evidence
#   "x"  = 6, # Exonic overlap on opposite strand
#   "o"  = 6, # Generic exonic overlap, no shared junctions
#   "i"  = 7, # Fully intronic, falls entirely within an intron of a reference transcript
#   "y"  = 8  # Contains reference within introns
#   )

# ── CATEGORISE MATCHES ──────────────────────────────────────
trmap_df_processed <- trmap_df_filt %>%
  mutate(
    cmp_category = case_when(
      cmp_class == "=" ~ "full_match",
      cmp_class == "k" & strand.sr == strand.lr ~ "same_loci" ,
      cmp_class %in% c("j", "c", "m", "n") ~ "same_loci",
      cmp_class %in% c("k", "x", "o", "i", "y") ~ "other_match",
      TRUE ~ "other_match"
    )
  )


# ── DEFINE TRANSCRIPT GROUPS ──────────────────────────────────────
lr_overlap <- unique(trmap_df_filt$transcript_id.lr)
sr_overlap <- unique(trmap_df_filt$transcript_id.sr)

all_lr <- all_transcripts %>%
          filter(num_samples_TPM_threshold.lr >= min_occurrence |
                  transcript_id.lr %in% lr_overlap) %>%
          pull(transcript_id.lr) %>%
          unique()
all_sr <- all_transcripts %>%
          filter(num_samples_TPM_threshold.sr >= min_occurrence |
                  transcript_id.sr %in% sr_overlap) %>%
          pull(transcript_id.sr) %>% unique()

lr_unique <- setdiff(all_lr, lr_overlap)
sr_unique <- setdiff(all_sr, sr_overlap)

# ── ASSES FULL MATCHES (cmp_class "=") ──────────────────────────────────────
full_matches <- trmap_df_processed %>%
  filter(cmp_category == "full_match") %>%
  mutate(type.lr = "transcript", 
        type.sr = "transcript")

# ── CALCULATE EXON METRICS ──────────────────────────────────────

# HELPERS

# Parse "start-end" string → numeric vector c(start, end)
parse_exon <- function(exon_str) {
  as.numeric(str_split(exon_str, "-")[[1]])
}

# Exon length from "start-end" string
exon_length <- function(exon_str) {
  coords <- parse_exon(exon_str)
  coords[2] - coords[1] + 1
}

# Extract terminal exon string from exon chain, strand-aware
# type = "first" → TSS-proximal exon
# type = "last"  → TES-proximal exon
# Exon chain is always in genomic (ascending) order:
#   + strand: first exon = chain[1],  last exon = chain[n]
#   - strand: first exon = chain[n],  last exon = chain[1]
terminal_exon <- function(exon_chain, strand, type) {
  exons <- str_split(exon_chain, ",")[[1]]
  if (strand == "+") {
    if (type == "first") exons[1] else tail(exons, 1)
  } else {
    if (type == "first") tail(exons, 1) else exons[1]
  }
}

# Junction position of a terminal exon (the internal splice site boundary)
# + strand first exon → end coordinate (right boundary = donor)
# + strand last exon  → start coordinate (left boundary = acceptor)
# - strand flipped
junction_pos <- function(exon_str, strand, type) {
  coords <- parse_exon(exon_str)
  if (strand == "+") {
    if (type == "first") coords[1] else coords[2]
  } else {
    if (type == "first") coords[2] else coords[1]
  }
}

# APPLY

full_matches_junc <- full_matches %>%
  rowwise() %>%
  mutate(
    exon_num.lr = str_count(exon_chain.lr, ",") + 1,
    exon_num.sr = str_count(exon_chain.sr, ",") + 1,
    
    first_exon.lr = terminal_exon(exon_chain.lr, strand.lr, "first"),
    first_exon.sr = terminal_exon(exon_chain.sr, strand.sr, "first"),
    last_exon.lr  = terminal_exon(exon_chain.lr, strand.lr, "last"),
    last_exon.sr  = terminal_exon(exon_chain.sr, strand.sr, "last"),
    
    first_exon_length.lr = exon_length(first_exon.lr),
    first_exon_length.sr = exon_length(first_exon.sr),
    last_exon_length.lr  = exon_length(last_exon.lr),
    last_exon_length.sr  = exon_length(last_exon.sr),
    
    lr_first_junc_pos = junction_pos(first_exon.lr, strand.lr, "first"),
    lr_last_junc_pos  = junction_pos(last_exon.lr,  strand.lr, "last"),
    sr_first_junc_pos = junction_pos(first_exon.sr, strand.sr, "first"),
    sr_last_junc_pos  = junction_pos(last_exon.sr,  strand.sr, "last")
  ) %>%
  ungroup()

# BUILD REFERENCE TERMINAL JUNCTIONS GRanges 
# reference is already defined in your script as lr_gtf_df with source != StringTie

ref_exons <- reference %>%
  filter(type == "exon") %>%
  makeGRangesFromDataFrame(
    seqnames.field     = "seqnames",
    start.field        = "start",
    end.field          = "end",
    strand.field       = "strand",
    keep.extra.columns = TRUE      # keeps gene_id, transcript_id
  )

# For each reference transcript extract the terminal junction positions
ref_terminal <- ref_exons %>%
  as.data.frame() %>%
  group_by(gene_id, transcript_id, seqnames, strand) %>%
  filter(exon_number == 1 | exon_number == max(exon_number)) %>%
  reframe(
    first_junc = if (unique(strand) == "+") start[exon_number   == 1]
    else                      end[exon_number   == 1],
    last_junc  = if (unique(strand) == "+") end[exon_number   != 1]
    else                         start[exon_number   != 1],
  )

# Collapse to gene level — one row per gene with vectors of annotated junctions
ref_junc_by_gene <- ref_terminal %>%
  group_by(gene_id, seqnames, strand) %>%
  summarise(
    first_juncs = list(unique(first_junc)),
    last_juncs  = list(unique(last_junc)),
  )

# ── RESOLVE TERMINAL EXONS ───────────────────────────────────────────────────

REF_MATCH_THRESHOLD <- 100  # nt — min LR/SR difference to trigger override
#      also used as max distance to ref to trust SR

closest_dist <- function(pos, ref_vec) {
  if (is.null(ref_vec) || length(ref_vec) == 0 || is.na(pos)) return(NA)
  min(abs(ref_vec - pos))
}

terminal_cmp <- full_matches_junc %>%
  select(
    transcript_id.lr, transcript_id.sr,
    gene_id.lr, strand.lr, seqnames.lr,
    first_exon.sr, last_exon.sr,
    first_exon.lr, last_exon.lr,
    first_exon_length.lr, first_exon_length.sr,
    last_exon_length.lr,  last_exon_length.sr,
    lr_first_junc_pos, lr_last_junc_pos,
    sr_first_junc_pos, sr_last_junc_pos
  ) %>%
  left_join(
    ref_junc_by_gene %>% dplyr::rename(gene_id.lr = gene_id),
    by = "gene_id.lr"
  ) %>%
  mutate(
    gene_in_reference    = !is.na(first_juncs),
    first_length_diff    = abs(first_exon_length.lr - first_exon_length.sr),
    last_length_diff     = abs(last_exon_length.lr  - last_exon_length.sr),
    
    first_exon_decision = case_when(
      lr_first_junc_pos  == sr_first_junc_pos       ~ "identical",
      first_length_diff  <= REF_MATCH_THRESHOLD    ~ "equivalent",
      # Large difference and gene is known — pick the end closer to reference
      gene_in_reference ~ {
        lr_d <- mapply(closest_dist, lr_first_junc_pos, first_juncs)
        sr_d <- mapply(closest_dist, sr_first_junc_pos, first_juncs)
        ifelse(lr_d <= sr_d, "use_lr", "sr_closer_to_ref")
      },
      # Large difference but novel gene — keep LR (longer is not necessarily better
      # without a reference to anchor to)
      TRUE                                                                   ~ "use_lr"
    ),
    
    last_exon_decision = case_when(
      lr_last_junc_pos  == sr_last_junc_pos       ~ "identical",
      last_length_diff  <= REF_MATCH_THRESHOLD    ~ "equivalent",
      gene_in_reference ~ {
        lr_d <- mapply(closest_dist, lr_last_junc_pos, last_juncs)
        sr_d <- mapply(closest_dist, sr_last_junc_pos, last_juncs)
        ifelse(lr_d <= REF_MATCH_THRESHOLD | sr_d <= REF_MATCH_THRESHOLD, 
              ifelse(lr_d <= sr_d, "use_lr", "sr_closer_to_ref"), 
              "use_lr")
      },
      TRUE                                                                   ~ "use_lr"
    ),
    
    model_decision = case_when(
      first_exon_decision == "equal"            &
        last_exon_decision  == "equal"                                       ~ "identical",
      first_exon_decision == "sr_closer_to_ref" &
        last_exon_decision  == "sr_closer_to_ref"                            ~ "use_sr",
      first_exon_decision == "sr_closer_to_ref" &
        last_exon_decision  != "sr_closer_to_ref"                            ~ "patch_5prime",
      first_exon_decision != "sr_closer_to_ref" &
        last_exon_decision  == "sr_closer_to_ref"                            ~ "patch_3prime",
      TRUE                                                                   ~ "use_lr"
    ),
    
    decision_confidence = case_when(
      gene_in_reference  & model_decision %in% c("use_lr", "use_sr",
                                                "identical")               ~ "high",
      gene_in_reference  & model_decision %in% c("patch_5prime",
                                                "patch_3prime")            ~ "medium",
      !gene_in_reference                                                      ~ "low",
      TRUE                                                                    ~ "medium"
    )
  )
terminal_cmp <- terminal_cmp %>%
  mutate(
    first_length_diff = first_exon_length.lr - first_exon_length.sr,
    last_length_diff  = last_exon_length.lr  - last_exon_length.sr)



# ── PATCH TERMINAL EXONS IN LR GTF DATAFRAME ────────────────────────────────
# For each full_match transcript, adjust exon coordinates according to
# model_decision before building the final GTF

patch_decisions <- terminal_cmp %>%
  select(transcript_id.lr, transcript_id.sr,
        model_decision,
        first_exon.lr, first_exon.sr,
        last_exon.lr,  last_exon.sr,
        strand.lr) %>%
  left_join(full_matches %>%
            select(transcript_id.lr, transcript_id.sr,
                  num_samples.sr,
                  num_samples_TPM_threshold.sr))

# Helper: parse "start-end" string → c(start, end)
parse_exon_str <- function(s) as.numeric(str_split(s, "-")[[1]])

lr_patched_gtf <- lr_gtf_df %>%
  right_join(patch_decisions) %>%
  rowwise() %>%
  mutate(
    # ── derive which exon is terminal at each end for this row ──────────
    is_first_exon = !is.na(model_decision) & type.lr == "exon" & {
      fe_coords <- parse_exon_str(first_exon.lr)
      start.lr == fe_coords[1] & end.lr == fe_coords[2]
    },
    is_last_exon = !is.na(model_decision) & type.lr == "exon" & {
      le_coords <- parse_exon_str(last_exon.lr)
      start.lr == le_coords[1] & end.lr == le_coords[2]
    },
    
    # ── apply patches ───────────────────────────────────────────────────
    start.lr = case_when(
      # patch_5prime: replace first exon start/end with SR first exon
      is_first_exon & model_decision == "patch_5prime" ~ {
        parse_exon_str(first_exon.sr)[1]
      },
      # patch_3prime: replace last exon with SR last exon
      is_last_exon  & model_decision == "patch_3prime" ~ {
        parse_exon_str(last_exon.sr)[1]
      },
      # use_sr: replace both terminal exons
      is_first_exon & model_decision == "use_sr" ~ parse_exon_str(first_exon.sr)[1],
      is_last_exon  & model_decision == "use_sr" ~ parse_exon_str(last_exon.sr)[1],
      TRUE ~ start.lr
    ),
    end.lr = case_when(
      is_first_exon & model_decision == "patch_5prime" ~ parse_exon_str(first_exon.sr)[2],
      is_last_exon  & model_decision == "patch_3prime" ~ parse_exon_str(last_exon.sr)[2],
      is_first_exon & model_decision == "use_sr"       ~ parse_exon_str(first_exon.sr)[2],
      is_last_exon  & model_decision == "use_sr"       ~ parse_exon_str(last_exon.sr)[2],
      TRUE ~ end.lr
    )
  ) %>%
  ungroup() %>%
  # also fix transcript-level start/end to match patched exons
  group_by(transcript_id.lr) %>%
  mutate(
    start.lr = ifelse(type.lr == "transcript", min(start.lr[type.lr == "exon"], na.rm = TRUE), start.lr),
    end.lr   = ifelse(type.lr == "transcript", max(end.lr[type.lr   == "exon"], na.rm = TRUE), end.lr)
  ) %>%
  ungroup() %>%
  dplyr::select(-is_first_exon, -is_last_exon,
        -first_exon.lr, -first_exon.sr,
        -last_exon.lr,  -last_exon.sr,
        -transcript_id.sr) %>%
  mutate(
    # NOTE: swapped to dplyr::if_else() — base ifelse() returns logical(0)
    # (instead of character(0)) when the input is zero-length, which breaks
    # the final bind_rows() if this piece ever has 0 rows.
    gene_id.lr   = dplyr::if_else(grepl("XLOC", gene_id.lr),
                          paste0(gene_id.lr, "_LR"), gene_id.lr),
    gene_name.lr = dplyr::if_else(grepl("XLOC", gene_name.lr),
                          paste0(gene_name.lr, "_LR"), gene_name.lr),
    transcript_id.lr    = paste0(transcript_id.lr, "_LR"),
    transcript_evidence.lr = dplyr::if_else(type.lr =="transcript", "LR_novel|SR_novel", NA_character_)) %>%
    rename_with(~ gsub("\\.lr$", "", .x),
              .cols = -any_of(c("num_samples_TPM_threshold.lr",
                                "num_samples.lr",
                                "num_samples_TPM_threshold.sr",
                                "num_samples.sr")))

# ── DEFINE lr_filtered AND sr_filtered ──────────────────────────────────────

same_loci_lr_ids <- trmap_df_processed %>%
  filter(cmp_category == "same_loci") %>%
  pull(transcript_id.lr) %>%
  unique()

same_loci_sr_ids <- trmap_df_processed %>%
  filter(cmp_category == "same_loci") %>%
  pull(transcript_id.sr) %>%
  unique()

lr_filtered <- lr_gtf_df %>%
  dplyr::filter(
    source.lr == "StringTie",
    transcript_id.lr %in% same_loci_lr_ids,
    ! transcript_id.lr %in% full_matches$transcript_id.lr,
    as.numeric(num_samples_TPM_threshold.lr) >= min_occurrence
  )

sr_filtered <- sr_gtf_df %>%
  dplyr::filter(
    source.sr == "StringTie",
    transcript_id.sr %in% same_loci_sr_ids,
    ! transcript_id.sr %in% full_matches$transcript_id.sr,
    as.numeric(num_samples_TPM_threshold.sr) >= min_occurrence
  )

# ── HARMONIZE GENE ANNOTATION FOR same_loci ──────────────────────────────────
# These are NOT merged — LR stays LR_novel, SR stays SR_novel
# Only gene_id and gene_name are unified so quantification works at gene level

gene_map <-trmap_df_processed %>%
  dplyr::filter(cmp_category %in% c("same_loci"), 
                !transcript_id.lr %in% full_matches$transcript_id.lr,
                !transcript_id.sr %in% full_matches$transcript_id.sr,
                grepl("XLOC", gene_name.lr)) %>%
  dplyr::select(transcript_id.lr, transcript_id.sr,
                gene_id.lr, gene_id.sr,
                gene_name.lr, gene_name.sr, cmp_class) %>%
  dplyr::distinct()

# For LR same_loci: adopt SR gene_id/gene_name if SR gene is known,
# otherwise keep LR gene with _LR suffix if XLOC
sr_same_loci_gtf <- sr_filtered %>%
  left_join(
    gene_map %>% select(transcript_id.sr, gene_id.lr, gene_name.lr) %>%
      distinct()
  ) %>%
  mutate(
    # NOTE: swapped to dplyr::if_else() — see comment above re: 0-row inputs
    gene_id.sr   = dplyr::if_else(!is.na(gene_id.lr),   paste0(gene_id.lr, "_LR"),
                          dplyr::if_else(grepl("XLOC", gene_id.sr),
                                paste0(gene_id.sr, "_SR"), gene_id.sr)),
    gene_name.sr = dplyr::if_else(!is.na(gene_name.lr),   paste0(gene_name.lr, "_LR"),
                          dplyr::if_else(grepl("XLOC", gene_id.sr),
                                paste0(gene_name.sr, "_SR"), gene_name.sr)),
    transcript_id.sr    = paste0(transcript_id.sr, "_SR"),
    transcript_evidence.sr = dplyr::if_else(type.sr == "transcript", transcript_evidence.sr, NA_character_),
    num_samples_TPM_threshold.lr = NA,
    num_samples.lr               = NA
  ) %>%
  select(-gene_id.lr, -gene_name.lr) %>%
  rename_with(~ gsub("\\.sr$", "", .x),
              .cols = -any_of(c("num_samples_TPM_threshold.lr",
                                "num_samples.lr",
                                "num_samples_TPM_threshold.sr",
                                "num_samples.sr")))

# For SR same_loci: adopt SR gene_id/gene_name (already has it),
# but also apply _SR suffix to XLOC gene IDs for consistency
lr_same_loci_gtf <- lr_filtered  %>%
  left_join(
    gene_map %>%
      select(transcript_id.lr, gene_id.lr, gene_name.lr) %>%
      distinct()) %>%
  mutate(
    # NOTE: swapped to dplyr::if_else() — see comment above re: 0-row inputs
    gene_id.lr   = dplyr::if_else(grepl("XLOC", gene_id.lr),
                          paste0(gene_id.lr, "_LR"), gene_id.lr),
    gene_name.lr = dplyr::if_else(grepl("XLOC", gene_name.lr),
                          paste0(gene_name.lr, "_LR"), gene_name.lr),
    transcript_id.lr    = paste0(transcript_id.lr, "_LR"),
    transcript_evidence.lr = dplyr::if_else(type.lr == "transcript", transcript_evidence.lr, NA_character_),
    num_samples_TPM_threshold.sr = NA,
    num_samples.sr               = NA
  ) %>%
  rename_with(~ gsub("\\.lr$", "", .x),
              .cols = -any_of(c("num_samples_TPM_threshold.sr",
                                "num_samples.sr",
                                "num_samples_TPM_threshold.lr",
                                "num_samples.lr")))

# ── LR-ONLY (no SR match at all) ─────────────────────────────────────────────

lr_rec_pass <- lr_gtf_df %>%
  filter(type.lr == "transcript", num_samples_TPM_threshold.lr >= min_occurrence) %>%
  pull(transcript_id.lr)

lr_only_gtf <- lr_gtf_df %>%
  filter(
    source.lr == "StringTie",
    ! transcript_id.lr %in% full_matches$transcript_id.lr, 
    ! transcript_id.lr %in% lr_filtered$transcript_id.lr,
    transcript_id.lr %in% lr_rec_pass
  ) %>%
  mutate(
    # NOTE: swapped to dplyr::if_else() — see comment above re: 0-row inputs
    gene_id.lr   = dplyr::if_else(grepl("XLOC", gene_id.lr),
                          paste0(gene_id.lr, "_LR"), gene_id.lr),
    gene_name.lr = dplyr::if_else(grepl("XLOC", gene_name.lr),
                          paste0(gene_name.lr, "_LR"), gene_name.lr),
    transcript_id.lr    = paste0(transcript_id.lr, "_LR"),
    transcript_evidence.lr = dplyr::if_else(type.lr == "transcript", transcript_evidence.lr, NA_character_),
    num_samples_TPM_threshold.sr = NA,
    num_samples.sr               =  NA
  ) %>%
  rename_with(~ gsub("\\.lr$", "", .x),
              .cols = -any_of(c("num_samples_TPM_threshold.lr",
                                "num_samples.lr",
                                "num_samples_TPM_threshold.sr",
                                "num_samples.sr")))

# ── SR-ONLY (no LR match at all) ─────────────────────────────────────────────

sr_rec_pass <- sr_gtf_df %>%
  filter(type.sr == "transcript", num_samples_TPM_threshold.sr >= min_occurrence) %>%
  pull(transcript_id.sr)

sr_only_gtf <- sr_gtf_df %>%
  filter(
    source.sr == "StringTie",
    ! transcript_id.sr %in% full_matches$transcript_id.sr, 
    ! transcript_id.sr %in% sr_filtered$transcript_id.sr, 
    transcript_id.sr %in% sr_rec_pass
  ) %>%
  mutate(
    # NOTE: swapped to dplyr::if_else() — this is the block that actually
    # crashed: with 0 rows, base ifelse() returned a logical(0) gene_id
    # column instead of a character(0) one, which vctrs::bind_rows()
    # refused to combine with the character gene_id columns from the
    # other pieces further down in the script.
    gene_id.sr   = dplyr::if_else(grepl("XLOC", gene_id.sr),
                          paste0(gene_id.sr, "_SR"), gene_id.sr),
    gene_name.sr = dplyr::if_else(grepl("XLOC", gene_name.sr),
                          paste0(gene_name.sr, "_SR"), gene_name.sr),
    transcript_id.sr    = paste0(transcript_id.sr, "_SR"),
    transcript_evidence.sr = dplyr::if_else(type.sr == "transcript", transcript_evidence.sr, NA_character_),
    num_samples_TPM_threshold.lr = NA,
    num_samples.lr               =NA
  ) %>%
  rename_with(~ gsub("\\.sr$", "", .x),
              .cols = -any_of(c("num_samples_TPM_threshold.sr",
                                "num_samples.sr",
                                "num_samples_TPM_threshold.lr",
                                "num_samples.lr")))

# ── REFERENCE — pass through unchanged ───────────────────────────────────────
reference_gtf <- reference %>%
  left_join(lr_gtf_df %>%
              dplyr::select(transcript_id.lr, gene_id.lr, type.lr, transcript_evidence.lr) %>% 
              filter(type.lr == "transcript") %>%
              distinct(), 
            by=join_by(transcript_id == transcript_id.lr, 
                      gene_id == gene_id.lr,
                      type == type.lr)) %>%
  mutate(
    transcript_evidence          = ifelse(type == "transcript", 
                                          paste(transcript_evidence.lr, transcript_evidence, sep = "|"), 
                                          NA),
    model_decision               = NA,
    num_samples_TPM_threshold.lr = NA,
    num_samples_TPM_threshold.sr = NA,
    num_samples.lr               = NA,
    num_samples.sr               = NA
  ) %>%
  select(-num_samples, -num_samples_TPM_threshold, -transcript_evidence.lr)

# ── BIND ALL PARTS ───────────────────────────────────────────────────────────
unified_gtf_df <- bind_rows(
  lr_patched_gtf %>% mutate(across(starts_with("num_"), as.numeric)),
  lr_same_loci_gtf %>% mutate(across(starts_with("num_"), as.numeric)),
  sr_same_loci_gtf %>% mutate(across(starts_with("num_"), as.numeric)),
  lr_only_gtf %>% mutate(across(starts_with("num_"), as.numeric)),
  sr_only_gtf %>% mutate(across(starts_with("num_"), as.numeric)),
  reference_gtf %>% mutate(across(starts_with("num_"), as.numeric))
) %>%
  mutate(
    across(c(num_samples_TPM_threshold.lr, num_samples_TPM_threshold.sr,
            num_samples.lr, num_samples.sr,
            transcript_evidence, model_decision),
          ~ ifelse(type == "transcript", .x, NA))
  ) %>%
  mutate(gene_biotype = ifelse(grepl("XLOC", gene_id), "StringTie", gene_biotype))

col_order <-  c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "phase",
                "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_biotype", "exon_number", 
                "class_code", "contained_in", "oId", "cmp_ref", "cmp_ref_gene","tss_id",  
                "transcript_evidence",  "model_decision", "num_samples.lr", "num_samples_TPM_threshold.lr", 
                "num_samples.sr", "num_samples_TPM_threshold.sr")

unified_gtf_df <- unified_gtf_df[,col_order]

unified_gtf_df_sorted <- sort_gtf(unified_gtf_df) 


# ── EXPORT TO GTF ─────────────────────────────────────────────────────────────
unified_gr <- makeGRangesFromDataFrame(
  unified_gtf_df_sorted,
  seqnames.field     = "seqnames",
  start.field        = "start",
  end.field          = "end",
  strand.field       = "strand",
  keep.extra.columns = TRUE
)

unified_gr_novel <- makeGRangesFromDataFrame(
  unified_gtf_df_sorted[unified_gtf_df_sorted$source == "StringTie",],
  seqnames.field     = "seqnames",
  start.field        = "start",
  end.field          = "end",
  strand.field       = "strand",
  keep.extra.columns = TRUE
)

output_gtf_path_novel <- file.path(qc_output, "unified_lr_sr_novel.gtf")
output_gtf_path <- file.path(qc_output, "unified_lr_sr.gtf")

## Write GTF files using unified helper function
write_gtf_with_header(unified_gr_novel, output_gtf_path_novel, gtf_header)
write_gtf_with_header(unified_gr, output_gtf_path, gtf_header)

# ── SUMMARY STATISTICS ───────────────────────────────────────────────────────

tx_only <- unified_gtf_df %>% filter(type == "transcript")

# ── Gene-level novelty ───────────────────────────────────────────────────────
# Novel loci = genes with XLOC gene_id
# Categorise each novel gene by the evidence supporting it:
#   LR_novel|SR_novel : gene has at least one fully matched transcript (both platforms)
#   LR_novel|SR_novel (mixed): gene has both LR_novel and SR_novel transcripts (no full match)
#   LR_novel only     : all transcripts under this gene are LR_novel
#   SR_novel only     : all transcripts under this gene are SR_novel

novel_gene_summary <- tx_only %>%
  filter(grepl("XLOC", gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    has_lr_sr_match = any(transcript_evidence == "LR_novel|SR_novel", na.rm = TRUE),
    has_lr          = any(transcript_evidence == "LR_novel",          na.rm = TRUE),
    has_sr          = any(transcript_evidence == "SR_novel",          na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    gene_novelty = case_when(
      has_lr_sr_match              ~ "LR_novel|SR_novel",
      has_lr & has_sr              ~ "LR_novel|SR_novel (mixed)",
      has_lr & !has_sr             ~ "LR_novel",
      !has_lr & has_sr             ~ "SR_novel"
    )
  )

gene_counts <- novel_gene_summary %>%
  count(gene_novelty, name = "n_genes")

# ── Transcript-level novelty ─────────────────────────────────────────────────
tx_counts <- tx_only %>%
  filter(grepl("TCONS", transcript_id)) %>%
  mutate(source_loci = ifelse(grepl("XLOC", gene_id), "novel", "annotated")) %>%
  group_by(source_loci) %>%
  count(transcript_evidence, name = "n_transcripts") %>%
  dplyr::rename(category = transcript_evidence)

# ── Log messages ─────────────────────────────────────────────────────────────
write_lines(c(
  "Unified GTF written\n",
  "  Total rows        : ", nrow(unified_gtf_df), "\n",
  "  Total transcripts : ", nrow(tx_only), "\n",
  "\n── Novel loci ──\n",
  "  Total novel loci            : ", nrow(novel_gene_summary), "\n",
  "  LR_novel|SR_novel           : ", sum(gene_counts$n_genes[gene_counts$gene_novelty == "LR_novel|SR_novel"],          na.rm = TRUE), "\n",
  "  LR_novel|SR_novel (mixed)   : ", sum(gene_counts$n_genes[gene_counts$gene_novelty == "LR_novel|SR_novel (mixed)"],  na.rm = TRUE), "\n",
  "  LR_novel only               : ", sum(gene_counts$n_genes[gene_counts$gene_novelty == "LR_novel"],                   na.rm = TRUE), "\n",
  "  SR_novel only               : ", sum(gene_counts$n_genes[gene_counts$gene_novelty == "SR_novel"],                   na.rm = TRUE), "\n",
  "\n── Novel transcripts ──\n",
  "  Total novel transcripts     : ", sum(tx_counts$n_transcripts), "\n",
  "  ...in annotated genes       : ", sum(tx_counts$n_transcripts[tx_counts$source_loci == "annotated"],  na.rm = TRUE), "\n",
  "  ...in novel loci            : ", sum(tx_counts$n_transcripts[tx_counts$source_loci == "novel"],  na.rm = TRUE), "\n",
  "  LR_novel|SR_novel           : ", sum(tx_counts$n_transcripts[tx_counts$category == "LR_novel|SR_novel"],  na.rm = TRUE), "\n",
  "  LR_novel                    : ", sum(tx_counts$n_transcripts[tx_counts$category == "LR_novel"],           na.rm = TRUE), "\n",
  "  SR_novel                    : ", sum(tx_counts$n_transcripts[tx_counts$category == "SR_novel"],           na.rm = TRUE)), "\n\n",
  sep = "",
  file =  file.path(qc_output, "unified_transcriptome.log")
)