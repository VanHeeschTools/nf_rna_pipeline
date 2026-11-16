#!/usr/bin/env python3
"""Merge the long-read and short-read transcriptomes into one annotation.

    lr_sr_combine.py <SR gtf> <LR gtf> <trmap.tab> <output dir> <min_occurrence> \
                     [end_tolerance] [threads]

Inputs
    SR gtf, LR gtf     the extended_reference.gtf of each pipeline: the reference plus that platform's
                       novel transcripts (source = StringTie)
    trmap.tab          trmap -T output: one line per overlapping LR-SR pair, with its class code
    output dir         folder for the output files
    min_occurrence     the number of samples in which a transcript (or pair) must be expressed
    end_tolerance      REF_MATCH_THRESHOLD of the R script, in nt (default 100, see step 3)
    threads            number of threads for polars and polars-bio (default: all cores, or the
                       POLARS_MAX_THREADS environment variable). In Nextflow: ${task.cpus}

Outputs
    unified_lr_sr.gtf              the reference plus all unified novel transcripts
    unified_lr_sr_novel.gtf        only the novel transcripts
    unified_transcriptome.log      counts per gene and transcript category

Requires: polars, polars-bio.
"""
import argparse
import os
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sr_gtf", help="short-read extended_reference.gtf")
    parser.add_argument("lr_gtf", help="long-read extended_reference.gtf")
    parser.add_argument("trmap_tab", help="trmap -T output")
    parser.add_argument("output_dir")
    parser.add_argument("min_occurrence", type=int)
    parser.add_argument("end_tolerance", nargs="?", type=int, default=100,
                        help="REF_MATCH_THRESHOLD of the R script, in nt (default 100)")
    parser.add_argument("threads", nargs="?", type=int, default=None,
                        help="threads for polars and polars-bio (default: all cores)")
    return parser


def parse_arguments(argv=None):
    return build_parser().parse_args(argv)


# The arguments are parsed before polars is imported: polars fixes its thread pool size at import.
THREADS = parse_arguments(sys.argv[1:]).threads
if THREADS:
    os.environ["POLARS_MAX_THREADS"] = str(THREADS)

import polars as pl          # noqa: E402  (after POLARS_MAX_THREADS is set)
import polars_bio as pb      # noqa: E402

if THREADS:
    pb.set_option("datafusion.execution.target_partitions", str(THREADS))
pb.set_option(pb.POLARS_BIO_COORDINATE_SYSTEM_ZERO_BASED, False)

GTF_ATTRIBUTES = [
    "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_biotype", "exon_number",
    "class_code", "cmp_ref", "cmp_ref_gene", "contained_in", "oId", "tss_id", "xloc",
    "num_samples", "num_samples_TPM_threshold", "transcript_evidence",
]
TRMAP_TAB_COLUMNS = ["lr_transcript", "class_code", "matched_percent", "sr_transcript", "matched_junctions"]
KEPT_ROUTES = ["full match", "same locus", "long-read only", "short-read only"]
OUTPUT_COLUMNS = [
    "chrom", "start", "end", "strand", "source", "type", "score", "phase",
    "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_biotype", "exon_number",
    "class_code", "contained_in", "oId", "cmp_ref", "cmp_ref_gene", "tss_id",
    "transcript_evidence", "model_decision",
    "num_samples.lr", "num_samples_TPM_threshold.lr", "num_samples.sr", "num_samples_TPM_threshold.sr",
]
FEATURE_TYPE_PRIORITY = {"gene": 1, "transcript": 2, "exon": 3, "CDS": 4, "UTR": 4, "start_codon": 4}
STRANDED_FEATURE_TYPES = ["exon", "CDS", "UTR", "start_codon"]
FIXED_GTF_COLUMNS = ["chrom", "source", "type", "start", "end", "score", "strand", "phase"]


def log(message: str) -> None:
    """Progress messages go to stderr, so they end up in the Nextflow .command.err."""
    print(message, file=sys.stderr, flush=True)


# =====================================================================================================
# Read gtf functions
# =====================================================================================================
def scan_gtf_table(path) -> pl.LazyFrame:
    """Lazily scan a GTF: 8 standard columns + one column per attribute (empty when absent)."""
    gtf = pb.scan_gtf(str(path), attr_fields=GTF_ATTRIBUTES)
    columns = gtf.collect_schema().names()
    missing = [name for name in FIXED_GTF_COLUMNS if name not in columns]
    if missing:
        raise KeyError(f"pb.scan_gtf() returned no column(s) {missing}; got {columns}")
    absent = [name for name in GTF_ATTRIBUTES if name not in columns]
    return gtf.with_columns(
        *[pl.lit(None, dtype=pl.String).alias(name) for name in absent],
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.when(pl.col("strand").is_in(["+", "-"])).then(pl.col("strand")).otherwise(pl.lit("*")).alias("strand"),
    )


def with_platform(novel_rows: pl.DataFrame, platform: str) -> pl.DataFrame:
    """Add the platform and make ids unique: TCONS_1 -> TCONS_1-LR, XLOC_1 -> XLOC_1-LR (ENSG unchanged)."""
    suffix = f"-{platform}"
    return novel_rows.with_columns(
        pl.lit(platform).alias("platform"),
        (pl.col("transcript_id") + suffix).alias("transcript_id"),
        pl.when(pl.col("gene_id").str.contains("XLOC")).then(pl.col("gene_id") + suffix)
          .otherwise(pl.col("gene_id")).alias("gene_id"),
    )


def transcript_summary(novel_rows: pl.DataFrame) -> pl.DataFrame:
    """One row per novel transcript, with what the later steps decide on.

    Columns: platform, transcript and gene id (suffixed), strand, the two sample-count attributes
    as written by filter_annotate (copied to the output in step 6), n (the number of samples above
    the TPM threshold, as a number) and exon_count.
    """
    exon_counts = novel_rows.filter(pl.col("type") == "exon").group_by("transcript_id").len("exon_count")
    return (
        novel_rows
        .filter(pl.col("type") == "transcript")
        .select(
            "platform", "transcript_id", "gene_id", "strand", "num_samples", "num_samples_TPM_threshold",
            # FIX bug 1: the sample count is a number; missing counts as 0
            pl.col("num_samples_TPM_threshold").cast(pl.Int64, strict=False).fill_null(0).alias("n"),
        )
        .join(exon_counts, on="transcript_id", how="left", maintain_order="left")
        .with_columns(pl.col("exon_count").fill_null(0))
    )


def terminal_exons(novel_rows: pl.DataFrame) -> pl.DataFrame:
    """Helper: per transcript, its 5' (first) and 3' (last) exon, TSS, TES and their lengths.

    Strand-aware: on the - strand the first exon is the right-most one, the TSS its end coordinate.
    """
    per_transcript = (
        novel_rows.filter(pl.col("type") == "exon").group_by("transcript_id")
        .agg(
            pl.col("strand").first(),
            pl.col("start").sort_by("start").first().alias("left_start"),
            pl.col("end").sort_by("start").first().alias("left_end"),
            pl.col("start").sort_by("start").last().alias("right_start"),
            pl.col("end").sort_by("start").last().alias("right_end"),
        )
    )
    plus = pl.col("strand") == "+"
    return per_transcript.select(
        "transcript_id",
        pl.when(plus).then(pl.col("left_start")).otherwise(pl.col("right_start")).alias("first_exon_start"),
        pl.when(plus).then(pl.col("left_end")).otherwise(pl.col("right_end")).alias("first_exon_end"),
        pl.when(plus).then(pl.col("right_start")).otherwise(pl.col("left_start")).alias("last_exon_start"),
        pl.when(plus).then(pl.col("right_end")).otherwise(pl.col("left_end")).alias("last_exon_end"),
        plus.alias("plus_strand"),
    ).with_columns(
        pl.when(pl.col("plus_strand")).then(pl.col("first_exon_start")).otherwise(pl.col("first_exon_end")).alias("tss"),
        pl.when(pl.col("plus_strand")).then(pl.col("last_exon_end")).otherwise(pl.col("last_exon_start")).alias("tes"),
        (pl.col("first_exon_end") - pl.col("first_exon_start") + 1).alias("first_exon_length"),
        (pl.col("last_exon_end") - pl.col("last_exon_start") + 1).alias("last_exon_length"),
    ).drop("plus_strand")


# =====================================================================================================
# Sorting and writing functions
# =====================================================================================================
def chromosome_sort_key() -> pl.Expr:
    """1..22, then X, Y, MT, then any other contig."""
    name = pl.col("chrom").str.replace(r"^chr", "")
    return (
        pl.when(name.str.contains(r"^\d+$")).then(name.cast(pl.Int64, strict=False))
          .when(name == "X").then(pl.lit(23))
          .when(name == "Y").then(pl.lit(24))
          .when(name.is_in(["M", "MT"])).then(pl.lit(25))
          .otherwise(pl.lit(1000))
    )


def sort_gtf(table: pl.DataFrame) -> pl.DataFrame:
    """Sort a GTF table gene by gene"""
    gene_rows = table.filter(pl.col("type") == "gene").select("chrom", "start", "end", "strand", "gene_id")
    genes_without_gene_row = (
        table
        .filter(pl.col("type") == "transcript")
        .join(gene_rows.select("gene_id"), on="gene_id", how="anti")
        .sort("start")
        .unique(subset="gene_id", keep="first", maintain_order=True)
        .select("chrom", "start", "end", "strand", "gene_id")
    )
    gene_rank = (
        pl.concat([gene_rows, genes_without_gene_row])
        .with_columns(chromosome_sort_key().alias("chromosome_key"))
        .sort(["chromosome_key", "chrom", "start", "end"], descending=[False, False, False, True])
        .unique(subset="gene_id", keep="first", maintain_order=True)
        .with_row_index("gene_rank")
        .select("gene_id", "gene_rank")
    )
    transcript_rank = (
        table
        .filter(pl.col("type") == "transcript")
        .select("gene_id", "transcript_id", "chrom", "start", "end")
        .unique()
        .with_columns(chromosome_sort_key().alias("chromosome_key"))
        .sort(["chromosome_key", "chrom", "start", "end", "gene_id", "transcript_id"],
              descending=[False, False, False, True, False, False])
        .unique(subset="transcript_id", keep="first", maintain_order=True)
        .with_columns((pl.int_range(pl.len()).over("gene_id") + 1).alias("transcript_rank"))
        .select("transcript_id", "transcript_rank")
    )
    alphabetical = sorted(table["type"].unique().to_list())
    rank_of_type = {feature: 4 + position for position, feature in enumerate(alphabetical, start=1)}
    rank_of_type.update(FEATURE_TYPE_PRIORITY)
    position_in_transcript = (
        pl.when(pl.col("type").is_in(STRANDED_FEATURE_TYPES) & (pl.col("strand") == "-"))
          .then(-pl.col("start")).otherwise(pl.col("start"))
    )
    return (
        table
        .join(gene_rank, on="gene_id", how="left", maintain_order="left")
        .join(transcript_rank, on="transcript_id", how="left", maintain_order="left")
        .with_columns(
            pl.col("type").replace_strict(rank_of_type, return_dtype=pl.Int64).alias("feature_rank"),
            position_in_transcript.alias("position_in_transcript"),
        )
        .sort(["gene_rank", "transcript_rank", "feature_rank", "position_in_transcript"],
              nulls_last=False, maintain_order=True)
        .drop("gene_rank", "transcript_rank", "feature_rank", "position_in_transcript")
    )


def write_gtf(table: pl.DataFrame, header_lines: list, path) -> None:
    """Write a GTF table: header lines, then 8 fixed columns + 'key "value";' attributes per row.
    """
    attribute_columns = [c for c in ("gene_id", "transcript_id") if c in table.columns]
    attribute_columns += [c for c in table.columns if c not in FIXED_GTF_COLUMNS + attribute_columns]
    attribute_text = pl.concat_str([
        pl.when(pl.col(column).is_not_null())
          .then(pl.lit(f'{column} "') + pl.col(column).cast(pl.String) + pl.lit('"; '))
          .otherwise(pl.lit(""))
        for column in attribute_columns
    ]).str.strip_chars_end(" ")
    lines = table.select(
        pl.concat_str(
            [
                pl.col("chrom"),
                pl.col("source").fill_null("."),
                pl.col("type"),
                pl.col("start").cast(pl.String),
                pl.col("end").cast(pl.String),
                pl.col("score").cast(pl.String).fill_null("."),
                pl.when(pl.col("strand") == "*").then(pl.lit(".")).otherwise(pl.col("strand")),
                pl.col("phase").cast(pl.String).fill_null("."),
                attribute_text,
            ],
            separator="\t",
        ).alias("line")
    )
    with open(path, "wb") as handle:
        handle.write(("\n".join(header_lines) + "\n").encode())
        lines.write_csv(handle, include_header=False, quote_style="never")


def connected_groups(links) -> list:
    """Split linked items into groups: items connected directly or through a chain of links share a group.

    links: (a, b) pairs, e.g. (long-read gene, short-read gene) for every kept LR-SR pair.
    Returns a list of sets, one per group. Example:

        [("A", "X"), ("B", "X"), ("C", "Y")]  ->  [{"A", "B", "X"}, {"C", "Y"}]

    A and B are not linked directly, but both are linked to X, so all three are one group.
    A simple "look up the partner" would miss this; that is what splits genes in the R script.

    Method (union-find): every item points to a parent; an item that points to itself is the
    representative of its group. Linking two items joins their groups by pointing one
    representative to the other. In the end, items with the same representative form one group.
    """
    parent = {}   # item -> parent item; new items start as their own group

    def representative(item):
        """Follow the parents up to the item that points to itself."""
        parent.setdefault(item, item)
        while parent[item] != item:
            parent[item] = parent[parent[item]]   # shortcut the path, so later lookups are faster
            item = parent[item]
        return item

    # join the groups of the two items of every link
    for a, b in links:
        representative_a, representative_b = representative(a), representative(b)
        if representative_a != representative_b:
            parent[representative_b] = representative_a

    # collect the items per representative
    groups = {}
    for item in list(parent):
        groups.setdefault(representative(item), set()).add(item)
    return list(groups.values())


# =====================================================================================================
# Functions
# =====================================================================================================
def read_inputs(sr_gtf_path, lr_gtf_path):
    """Read both extended_reference GTFs.

    Returns
        sr_reference           every non-StringTie row of the SR GTF: the reference that goes
                               into the output (both GTFs contain the same reference)
        lr_reference_evidence  per reference transcript, its LR evidence label (LR_full_match, …),
                               to combine with the SR label in step 5
        novel_rows             the StringTie rows (transcripts and exons) of both platforms in one
                               table, with a platform column and -LR / -SR suffixed ids
    """
    is_novel = pl.col("source") == "StringTie"
    is_transcript = pl.col("type") == "transcript"
    # every query gets its own fresh scan (polars-bio scans should not be shared between queries)
    sr_reference = scan_gtf_table(sr_gtf_path).filter(~is_novel).collect()
    lr_reference_evidence = (
        scan_gtf_table(lr_gtf_path).filter(~is_novel & is_transcript)
        .select("transcript_id", "gene_id", pl.col("transcript_evidence").alias("lr_evidence"))
        .collect()
    )
    novel_rows = pl.concat(
        [with_platform(scan_gtf_table(lr_gtf_path).filter(is_novel).collect(), "LR"),
         with_platform(scan_gtf_table(sr_gtf_path).filter(is_novel).collect(), "SR")],
        how="diagonal_relaxed",
    )
    return sr_reference, lr_reference_evidence, novel_rows


def check_inputs(sr_reference: pl.DataFrame, transcripts: pl.DataFrame) -> None:
    """Warn (stderr) about inputs that make the output look wrong without an error."""
    if sr_reference.filter(pl.col("type") == "transcript").is_empty():
        log("WARNING: the SR GTF has no reference transcripts (source != StringTie): is it the "
            "extended_reference.gtf and not the novel_transcripts.gtf?")
    for platform, total, missing in (
        transcripts.group_by("platform")
        .agg(pl.len().alias("total"), pl.col("num_samples_TPM_threshold").is_null().sum().alias("missing"))
        .iter_rows()
    ):
        if missing:
            log(f"WARNING: {missing} of {total} {platform} transcripts have no num_samples_TPM_threshold; "
                f"they count as n = 0 (filter_annotate only writes it with min_tpm > 0)")


def read_pairs(trmap_tab_path, transcripts: pl.DataFrame, min_occurrence: int) -> pl.DataFrame:
    """Every trmap pair, with its relation and whether it is recurrent.

    trmap -T writes one line per overlapping pair: LR transcript (query), class code, matched %,
    SR transcript (reference), matched junctions. The class code describes the LR transcript
    relative to the SR one, with the gffcompare codes. Per pair:

        relation   "full match"  class =                       the same intron chain
                   "same locus"  class j, c, m, n, or k on     related isoforms of one gene
                                 the same strand
                   "other"       anything else (e.g. x)        they only overlap
        recurrent  full match: n_LR + n_SR >= min_occurrence   one transcript seen twice
                   otherwise:  n_LR >= min_occurrence and      two transcripts, each must be
                               n_SR >= min_occurrence          recurrent on its own

    Pairs whose transcripts are not in the novel table are dropped (inner joins), and so are pairs
    with a single-exon LR transcript (those are dropped in step 2 anyway). matched % and matched
    junctions are only used to choose between several full-match partners in step 2.
    """
    def pair_side(platform: str, name: str) -> pl.DataFrame:
        return transcripts.filter(pl.col("platform") == platform).select(
            pl.col("transcript_id").alias(name), pl.col("gene_id").alias(f"{name}_gene"),
            pl.col("strand").alias(f"{name}_strand"), pl.col("n").alias(f"n_{name}"),
            pl.col("exon_count").alias(f"{name}_exons"),
        )

    raw_pairs = pl.read_csv(trmap_tab_path, separator="\t", has_header=False, quote_char=None,
                            infer_schema=False, new_columns=TRMAP_TAB_COLUMNS)
    same_locus_code = pl.col("class_code").is_in(["j", "c", "m", "n"]) | (
        (pl.col("class_code") == "k") & (pl.col("lr_strand") == pl.col("sr_strand")))
    return (
        raw_pairs
        .select(
            (pl.col("lr_transcript") + "-LR").alias("lr"),
            (pl.col("sr_transcript") + "-SR").alias("sr"),
            "class_code",
            # the first number of each score ("85.0", "3/4", …); only used to break ties
            pl.col("matched_percent").str.extract(r"([\d.]+)", 1).cast(pl.Float64),
            pl.col("matched_junctions").str.extract(r"(\d+)", 1).cast(pl.Int64),
        )
        .unique(subset=["lr", "sr"], keep="first", maintain_order=True)
        .join(pair_side("LR", "lr"), on="lr", how="inner", maintain_order="left")
        .join(pair_side("SR", "sr"), on="sr", how="inner", maintain_order="left")
        .filter(pl.col("lr_exons") >= 2)          # no long-read single-exon transcripts
        .with_columns(
            pl.when(pl.col("class_code") == "=").then(pl.lit("full match"))
              .when(same_locus_code).then(pl.lit("same locus"))
              .otherwise(pl.lit("other")).alias("relation")
        )
        .with_columns(
            # full match: one transcript seen twice, so the samples are pooled
            pl.when(pl.col("relation") == "full match").then(pl.col("n_lr") + pl.col("n_sr") >= min_occurrence)
              .otherwise((pl.col("n_lr") >= min_occurrence) & (pl.col("n_sr") >= min_occurrence))
              .alias("recurrent")
        )
    )


def decide_routes(pairs: pl.DataFrame, transcripts: pl.DataFrame, min_occurrence: int):
    """Give every novel transcript exactly one route.

    Only recurrent pairs count. The rules, in order; the first that applies decides:

        dropped          LR transcript with a single exon
        full match       LR side of a recurrent full match; with several = partners, the best one
                         (highest matched %, then most junctions) is its partner (FIX bug 5)
        merged           SR side of a recurrent full match: it becomes part of that LR model
        same locus       in a recurrent same-locus pair: kept as it is (its gene: step 4)
        long-read only / short-read only
                         no recurrent full-match or same-locus pair, but n >= min_occurrence
        dropped          everything else (too few samples, and no recurrent partner)

    Why single-exon LR transcripts are still here: the long-read filter_annotate only removes
    single-exon transcripts without a reference gene, the short-read one removes all of them.

    Returns
        routes            the transcripts table with a route column
        best_full_match   one row per full-match model: the LR transcript and its chosen SR partner
        recurrent_pairs   the recurrent pairs, for the gene links of step 4
    """
    recurrent_pairs = pairs.filter(pl.col("recurrent"))
    full_match_pairs = recurrent_pairs.filter(pl.col("relation") == "full match")
    same_locus_pairs = recurrent_pairs.filter(pl.col("relation") == "same locus")

    # FIX bug 5: one partner per long-read model, the best match
    best_full_match = (
        full_match_pairs
        .sort(["lr", "matched_percent", "matched_junctions", "sr"], descending=[False, True, True, False], nulls_last=True)
        .unique(subset="lr", keep="first", maintain_order=True)
    )

    # the transcripts in each kind of recurrent pair
    model_ids = full_match_pairs["lr"].unique().to_list()
    merged_ids = full_match_pairs["sr"].unique().to_list()
    same_locus_ids = pl.concat([same_locus_pairs["lr"], same_locus_pairs["sr"]]).unique().to_list()

    # the rules of the docstring, in the same order
    is_long_read = pl.col("platform") == "LR"
    route = (
        pl.when(is_long_read & (pl.col("exon_count") < 2)).then(pl.lit("dropped"))
          .when(pl.col("transcript_id").is_in(model_ids)).then(pl.lit("full match"))
          .when(pl.col("transcript_id").is_in(merged_ids)).then(pl.lit("merged"))
          .when(pl.col("transcript_id").is_in(same_locus_ids)).then(pl.lit("same locus"))
          .when(pl.col("n") >= min_occurrence)
          .then(pl.when(is_long_read).then(pl.lit("long-read only")).otherwise(pl.lit("short-read only")))
          .otherwise(pl.lit("dropped"))
    )
    routes = transcripts.with_columns(route.alias("route"))
    return routes, best_full_match, recurrent_pairs


def decide_ends(best_full_match, novel_rows, sr_reference, tolerance: int) -> pl.DataFrame:
    """For every full match, decide which 5' and 3' end the model gets.

    The model is the LR transcript: long reads usually span the whole molecule, and the intron
    chain is the same anyway. Only the outer ends (TSS, TES) can differ. Per end, with tolerance
    t (--end-tolerance, REF_MATCH_THRESHOLD in R), the first rule that applies decides:

        1  the LR and SR end are at the same position             identical   (keep LR)
        2  the terminal exons differ by at most t in length       equivalent  (keep LR)
        3  new gene: no annotated TSS / TES to compare with       use_lr
        4  3' end only: neither end is within t of an annotated   use_lr
           TES of the gene
        5  otherwise the end closest to an annotated TSS / TES    use_lr or use_sr
           of the same gene wins (a tie keeps LR)

    The two ends give the model decision: use_sr (both ends SR), patch_5prime or patch_3prime
    (one end SR), identical (both identical), or use_lr.

    The annotated TSS / TES come from the reference transcript rows (FIX bug 3); a gene counts as
    annotated when it has reference transcripts (FIX bug 6); 'identical' can occur (FIX bug 7).

    Returns one row per model (lr = its transcript id): the terminal exons of both transcripts,
    the distances to the nearest annotated end, decision_5prime, decision_3prime, model_decision.
    """
    # the LR and SR terminal exons of every full match side by side (columns ..._lr, ..._sr)
    ends = terminal_exons(novel_rows)
    end_comparison = (
        best_full_match.select("lr", "sr", "lr_gene")
        .join(ends.rename(lambda c: c if c == "transcript_id" else f"{c}_lr"),
              left_on="lr", right_on="transcript_id", how="left", maintain_order="left")
        .join(ends.rename(lambda c: c if c == "transcript_id" else f"{c}_sr"),
              left_on="sr", right_on="transcript_id", how="left", maintain_order="left")
    )
    # FIX bug 3: annotated TSS / TES from the reference transcript rows, not max(exon_number) on text
    plus = pl.col("strand") == "+"
    annotated_ends = sr_reference.filter(pl.col("type") == "transcript").select(
        "gene_id",
        pl.when(plus).then(pl.col("start")).otherwise(pl.col("end")).alias("annotated_tss"),
        pl.when(plus).then(pl.col("end")).otherwise(pl.col("start")).alias("annotated_tes"),
    )
    distances = (
        end_comparison.join(annotated_ends, left_on="lr_gene", right_on="gene_id", how="inner")
        .group_by("lr")
        .agg(
            (pl.col("annotated_tss") - pl.col("tss_lr")).abs().min().alias("tss_distance_lr"),
            (pl.col("annotated_tss") - pl.col("tss_sr")).abs().min().alias("tss_distance_sr"),
            (pl.col("annotated_tes") - pl.col("tes_lr")).abs().min().alias("tes_distance_lr"),
            (pl.col("annotated_tes") - pl.col("tes_sr")).abs().min().alias("tes_distance_sr"),
        )
    )
    end_comparison = end_comparison.join(distances, on="lr", how="left", maintain_order="left")

    def end_decision(position: str, exon: str, distance: str, three_prime: bool) -> pl.Expr:
        """The rules 1-5 of the docstring for one end (rule 4 only for the 3' end)."""
        same_position = pl.col(f"{position}_lr") == pl.col(f"{position}_sr")
        length_difference = (pl.col(f"{exon}_length_lr") - pl.col(f"{exon}_length_sr")).abs()
        lr_distance, sr_distance = pl.col(f"{distance}_lr"), pl.col(f"{distance}_sr")
        gene_is_annotated = lr_distance.is_not_null()   # FIX bug 6
        decision = (
            pl.when(same_position).then(pl.lit("identical"))
              .when(length_difference <= tolerance).then(pl.lit("equivalent"))
              .when(~gene_is_annotated).then(pl.lit("use_lr"))
        )
        if three_prime:
            decision = decision.when((lr_distance > tolerance) & (sr_distance > tolerance)).then(pl.lit("use_lr"))
        return decision.when(lr_distance <= sr_distance).then(pl.lit("use_lr")).otherwise(pl.lit("use_sr"))

    return end_comparison.with_columns(
        end_decision("tss", "first_exon", "tss_distance", three_prime=False).alias("decision_5prime"),
        end_decision("tes", "last_exon", "tes_distance", three_prime=True).alias("decision_3prime"),
    ).with_columns(
        pl.when((pl.col("decision_5prime") == "use_sr") & (pl.col("decision_3prime") == "use_sr")).then(pl.lit("use_sr"))
          .when(pl.col("decision_5prime") == "use_sr").then(pl.lit("patch_5prime"))
          .when(pl.col("decision_3prime") == "use_sr").then(pl.lit("patch_3prime"))
          .when((pl.col("decision_5prime") == "identical") & (pl.col("decision_3prime") == "identical")).then(pl.lit("identical"))
          .otherwise(pl.lit("use_lr"))
          .alias("model_decision"),
    )


def patch_models(decisions: pl.DataFrame, novel_rows: pl.DataFrame) -> pl.DataFrame:
    """The GTF rows (transcript + exons) of every full-match model.

    Each model is its LR transcript. Where decide_ends() chose the SR end, the LR first exon
    (patch_5prime, use_sr) and / or last exon (patch_3prime, use_sr) gets the coordinates of the
    SR one, and the transcript row is stretched to the new exon span. The rows keep the LR ids and
    get a model_decision column.
    """
    patch_info = decisions.select(
        pl.col("lr").alias("transcript_id"), "model_decision",
        "first_exon_start_lr", "first_exon_end_lr", "last_exon_start_lr", "last_exon_end_lr",
        "first_exon_start_sr", "first_exon_end_sr", "last_exon_start_sr", "last_exon_end_sr",
    )
    is_exon = pl.col("type") == "exon"
    is_first_exon = is_exon & (pl.col("start") == pl.col("first_exon_start_lr")) & (pl.col("end") == pl.col("first_exon_end_lr"))
    is_last_exon = is_exon & (pl.col("start") == pl.col("last_exon_start_lr")) & (pl.col("end") == pl.col("last_exon_end_lr"))
    patch_first = is_first_exon & pl.col("model_decision").is_in(["patch_5prime", "use_sr"])
    patch_last = is_last_exon & pl.col("model_decision").is_in(["patch_3prime", "use_sr"])
    return (
        novel_rows
        .join(patch_info, on="transcript_id", how="inner", maintain_order="left")
        .with_columns(
            pl.when(patch_first).then(pl.col("first_exon_start_sr"))
              .when(patch_last).then(pl.col("last_exon_start_sr"))
              .otherwise(pl.col("start")).alias("start"),
            pl.when(patch_first).then(pl.col("first_exon_end_sr"))
              .when(patch_last).then(pl.col("last_exon_end_sr"))
              .otherwise(pl.col("end")).alias("end"),
        )
        .with_columns(   # the transcript row spans its (possibly patched) exons
            pl.when(pl.col("type") == "transcript").then(pl.col("start").filter(is_exon).min().over("transcript_id"))
              .otherwise(pl.col("start")).alias("start"),
            pl.when(pl.col("type") == "transcript").then(pl.col("end").filter(is_exon).max().over("transcript_id"))
              .otherwise(pl.col("end")).alias("end"),
        )
        .drop([c for c in patch_info.columns if c not in ("transcript_id", "model_decision")])
    )


def group_genes(recurrent_pairs: pl.DataFrame) -> pl.DataFrame:
    """Give every new gene found by both platforms one gene id.

    Returns one row per linked gene: gene_id, unified_gene_id and genes_in_group (the group size;
    large groups can point to read-through transcripts). Genes without a kept pair are not listed
    and keep their own id.
    """
    # the links: the genes of both transcripts of every kept pair, when both genes are new
    is_new_gene = pl.col("lr_gene").str.contains("XLOC") & pl.col("sr_gene").str.contains("XLOC")
    gene_links = (
        recurrent_pairs
        .filter(pl.col("relation").is_in(["full match", "same locus"]) & is_new_gene)
        .select("lr_gene", "sr_gene").unique()
    )
    rows = []
    for group in connected_groups(gene_links.iter_rows()):
        # every link has a long-read gene, so every group has one; the fallback is only a safeguard
        long_read_genes = sorted(gene for gene in group if gene.endswith("-LR"))
        group_id = long_read_genes[0] if long_read_genes else sorted(group)[0]
        rows += [(gene, group_id, len(group)) for gene in sorted(group)]
    return pl.DataFrame(rows, schema={"gene_id": pl.String, "unified_gene_id": pl.String, "genes_in_group": pl.Int64},
                        orient="row")


def reference_with_evidence(sr_reference: pl.DataFrame, lr_reference_evidence: pl.DataFrame) -> pl.DataFrame:
    """The reference rows of the output, with the evidence of both platforms.

    filter_annotate labelled every reference transcript with its support on that platform
    (e.g. LR_full_match, SR_no_evidence). Here the two labels are joined as "LR label|SR label",
    e.g. LR_full_match|SR_partial_match. A missing label is left out instead of becoming
    "NA|…". Exon rows get no label.
    """
    lr_evidence = lr_reference_evidence.unique(subset=["transcript_id", "gene_id"], keep="first")
    return (
        sr_reference
        .join(lr_evidence, on=["transcript_id", "gene_id"], how="left", maintain_order="left")
        .with_columns(
            pl.when(pl.col("type") == "transcript")
              .then(pl.concat_str([pl.col("lr_evidence"), pl.col("transcript_evidence")], separator="|", ignore_nulls=True))
              .alias("transcript_evidence")
        )
        .with_columns(pl.when(pl.col("transcript_evidence") != "").then(pl.col("transcript_evidence")).alias("transcript_evidence"))
        .drop("lr_evidence")
    )


def assemble(routes, best_full_match, models, gene_groups, novel_rows, reference_rows, transcripts, labels):
    """Put the output together as one sorted GTF table.

    Three parts are stacked:
        the full-match models                     
        the other kept transcripts (routes same locus,     
        long-read only, short-read only)
        the reference rows
    Merged and dropped transcripts are left out. Every novel transcript then gets:
        - its gene id: the group id for linked new genes, otherwise its own (suffixed) id
        - gene_name: the gene id for new genes; annotated genes keep their name
        - gene_biotype "stringtie" for new genes
        - num_samples / num_samples_TPM_threshold of both platforms (.lr / .sr): its own counts,
          and for a model also those of its SR partner
    Only OUTPUT_COLUMNS are kept, and the table is sorted gene by gene (sort_gtf).
    """
    lr_novel, sr_novel, both_novel = labels
    counts = transcripts.select(
        "transcript_id", "platform",
        pl.col("num_samples").alias("own_num_samples"),
        pl.col("num_samples_TPM_threshold").alias("own_tpm_samples"),
    )
    partner_counts = best_full_match.select(pl.col("lr").alias("transcript_id"), "sr").join(
        counts.select(pl.col("transcript_id").alias("sr"),
                      pl.col("own_num_samples").alias("partner_num_samples"),
                      pl.col("own_tpm_samples").alias("partner_tpm_samples")),
        on="sr", how="left",
    ).drop("sr")

    kept = routes.filter(pl.col("route").is_in(KEPT_ROUTES)).select("transcript_id", "route")
    model_ids = models["transcript_id"].unique().to_list()
    novel_output = pl.concat(
        [models, novel_rows.filter(~pl.col("transcript_id").is_in(model_ids))], how="diagonal_relaxed",
    ).join(kept, on="transcript_id", how="inner", maintain_order="left")

    is_transcript_row = pl.col("type") == "transcript"
    from_long_read = pl.col("platform") == "LR"

    def on_transcript_rows(value: pl.Expr) -> pl.Expr:
        return pl.when(is_transcript_row).then(value)

    novel_output = (
        novel_output
        .join(counts.drop("platform"), on="transcript_id", how="left", maintain_order="left")
        .join(partner_counts, on="transcript_id", how="left", maintain_order="left")
        .join(gene_groups.select("gene_id", "unified_gene_id"), on="gene_id", how="left", maintain_order="left")
        .with_columns(pl.coalesce("unified_gene_id", "gene_id").alias("gene_id"))
        .with_columns(
            pl.when(pl.col("gene_name").is_null() | pl.col("gene_name").str.contains("XLOC"))
              .then(pl.col("gene_id")).otherwise(pl.col("gene_name")).alias("gene_name"),
            pl.when(pl.col("gene_id").str.contains("XLOC")).then(pl.lit("stringtie"))
              .otherwise(pl.col("gene_biotype")).alias("gene_biotype"),
            on_transcript_rows(
                pl.when(pl.col("route") == "full match").then(pl.lit(both_novel))
                  .when(from_long_read).then(pl.lit(lr_novel))
                  .otherwise(pl.lit(sr_novel))
            ).alias("transcript_evidence"),
            on_transcript_rows(pl.col("model_decision")).alias("model_decision"),
            on_transcript_rows(pl.when(from_long_read).then(pl.col("own_num_samples"))).alias("num_samples.lr"),
            on_transcript_rows(pl.when(from_long_read).then(pl.col("own_tpm_samples"))).alias("num_samples_TPM_threshold.lr"),
            on_transcript_rows(pl.when(from_long_read).then(pl.col("partner_num_samples"))
                               .otherwise(pl.col("own_num_samples"))).alias("num_samples.sr"),
            on_transcript_rows(pl.when(from_long_read).then(pl.col("partner_tpm_samples"))
                               .otherwise(pl.col("own_tpm_samples"))).alias("num_samples_TPM_threshold.sr"),
        )
    )
    stacked = pl.concat([novel_output, reference_rows], how="diagonal_relaxed")
    return sort_gtf(stacked.select([c for c in OUTPUT_COLUMNS if c in stacked.columns]))


def build_log(unified_gtf, labels, args) -> list:
    """The log: new genes and novel transcripts per evidence category.

    A new gene's category: LR-novel|SR-novel when it has a full-match model, "… (mixed)" when it
    has transcripts of both platforms but no model, otherwise LR-novel or SR-novel. The lines that
    start with ## are the header of both GTF files.
    """
    lr_novel, sr_novel, both_novel = labels
    novel_transcripts = unified_gtf.filter((pl.col("type") == "transcript") & (pl.col("source") == "StringTie"))
    gene_categories = (
        novel_transcripts.filter(pl.col("gene_id").str.contains("XLOC"))
        .group_by("gene_id")
        .agg(
            (pl.col("transcript_evidence") == both_novel).any().alias("has_full_match"),
            (pl.col("transcript_evidence") == lr_novel).any().alias("has_lr"),
            (pl.col("transcript_evidence") == sr_novel).any().alias("has_sr"),
        )
        .with_columns(
            pl.when(pl.col("has_full_match")).then(pl.lit(both_novel))
              .when(pl.col("has_lr") & pl.col("has_sr")).then(pl.lit(f"{both_novel} (mixed)"))
              .when(pl.col("has_lr")).then(pl.lit(lr_novel))
              .otherwise(pl.lit(sr_novel)).alias("gene_category")
        )
    )

    def count_of(table: pl.DataFrame, column: str, value: str) -> int:
        return table.filter(pl.col(column) == value).height

    in_annotated_genes = novel_transcripts.filter(~pl.col("gene_id").str.contains("XLOC")).height
    mixed = f"{both_novel} (mixed)"
    return [
        "##UNIFIED LR/SR GTF (Python/Polars port, bugs fixed)",
        f"##SR GTF={args.sr_gtf}",
        f"##LR GTF={args.lr_gtf}",
        f"##trmap file={args.trmap_tab}",
        f"##min_occurrence={args.min_occurrence}",
        "",
        "Unified GTF written",
        f"  Total rows        : {unified_gtf.height}",
        f"  Total transcripts : {unified_gtf.filter(pl.col('type') == 'transcript').height}",
        "",
        "── Novel loci ──",
        f"  Total novel loci            : {gene_categories.height}",
        f"  {both_novel:<27} : {count_of(gene_categories, 'gene_category', both_novel)}",
        f"  {mixed:<27} : {count_of(gene_categories, 'gene_category', mixed)}",
        f"  {lr_novel:<27} : {count_of(gene_categories, 'gene_category', lr_novel)}",
        f"  {sr_novel:<27} : {count_of(gene_categories, 'gene_category', sr_novel)}",
        "",
        "── Novel transcripts ──",
        f"  Total novel transcripts     : {novel_transcripts.height}",
        f"  ...in annotated genes       : {in_annotated_genes}",
        f"  ...in novel loci            : {novel_transcripts.height - in_annotated_genes}",
        f"  {both_novel:<27} : {count_of(novel_transcripts, 'transcript_evidence', both_novel)}",
        f"  {lr_novel:<27} : {count_of(novel_transcripts, 'transcript_evidence', lr_novel)}",
        f"  {sr_novel:<27} : {count_of(novel_transcripts, 'transcript_evidence', sr_novel)}",
    ]


# =====================================================================================================
# Main
# =====================================================================================================
def main(argv=None) -> int:
    args = parse_arguments(argv)
    labels = ("LR-novel", "SR-novel", "LR-novel|SR-novel")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log(f"polars uses {pl.thread_pool_size()} threads")
    log("step 0: reading both GTFs")
    sr_reference, lr_reference_evidence, novel_rows = read_inputs(args.sr_gtf, args.lr_gtf)
    transcripts = transcript_summary(novel_rows)
    check_inputs(sr_reference, transcripts)

    log("step 1: trmap pairs")
    pairs = read_pairs(args.trmap_tab, transcripts, args.min_occurrence)
    log("step 2: routes")
    routes, best_full_match, recurrent_pairs = decide_routes(pairs, transcripts, args.min_occurrence)
    log("step 3: full-match models")
    end_decisions = decide_ends(best_full_match, novel_rows, sr_reference, args.end_tolerance)
    models = patch_models(end_decisions, novel_rows)
    log("step 4: gene groups")
    gene_groups = group_genes(recurrent_pairs)
    log("step 5: reference evidence")
    reference_rows = reference_with_evidence(sr_reference, lr_reference_evidence)
    log("step 6: assemble and write")
    unified_gtf = assemble(routes, best_full_match, models, gene_groups, novel_rows, reference_rows, transcripts, labels)
    unified_novel_gtf = unified_gtf.filter(pl.col("source") == "StringTie")

    duplicates = unified_gtf.filter(pl.col("type") == "transcript").group_by("transcript_id").len().filter(pl.col("len") > 1)
    if not duplicates.is_empty():
        raise RuntimeError(f"{duplicates.height} transcript ids occur more than once in the unified GTF")

    log_lines = build_log(unified_gtf, labels, args)
    gtf_header = [line for line in log_lines if line.startswith("##")]
    write_gtf(unified_gtf, gtf_header, output_dir / "unified_lr_sr.gtf")
    write_gtf(unified_novel_gtf, gtf_header, output_dir / "unified_lr_sr_novel.gtf")
    (output_dir / "unified_transcriptome.log").write_text("\n".join(log_lines) + "\n")

    log(f"done: {unified_novel_gtf.filter(pl.col('type') == 'transcript').height} novel transcripts in the unified GTF")
    return 0


if __name__ == "__main__":
    sys.exit(main())
