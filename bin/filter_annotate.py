#!/usr/bin/env python3
"""Filter and annotate assembled transcripts, and add them to the reference annotation.

    filter_annotate.py <reference.gtf> <combined.gtf> <combined.tracking> \
                        <min_occurrence> <min_tpm> <output_prefix> [platform] [threads] [refseq_basename]

Inputs
    reference.gtf        the reference annotation (Ensembl or GENCODE)
    combined.gtf         gffcompare -i output: every distinct transcript of all samples, with a class code
    combined.tracking    gffcompare tracking file: the TPM of every transcript in every sample
    min_occurrence       minimum number of samples that must express a transcript
    min_tpm              TPM at which a sample counts as expressing a transcript (0: use num_samples)
    output_prefix        prefix of all output files
    platform             LR (default): labels LR_…, single-exon filter only for transcripts without a
                         reference gene (long-read script). SR: labels SR_…, every single-exon transcript
                         is removed (short-read script).
    threads              number of threads for polars and polars-bio (default: all cores, or the
                         POLARS_MAX_THREADS environment variable). In Nextflow: ${task.cpus}
    refseq_basename      optional: reads <basename>.xr.gff and <basename>.nr.gff and flags novel
                         transcripts that overlap them (step 8); nothing is removed. Needs platform and
                         threads before it.

Outputs
    <prefix>.extended_reference.gtf   the reference plus the added transcripts
    <prefix>.novel_transcripts.gtf    the added transcripts only (not written when nothing is added)
    <prefix>.tsv                      TPM per sample of the added transcripts
    <prefix>.log                      counts per step

Steps
    0  read; keep assembled transcripts on the main chromosomes, with a strand
    1  label every reference transcript with its evidence (full / partial / other / no evidence)
    2  remove class codes = c j m n e r s (structure already in the reference)
    3  attach the TPM per sample from the tracking file, and assign gene ids (with the correction)
    4  remove single-exon transcripts
    5  remove transcripts expressed in too few samples (with the correction)
    6  assign a gene biotype, and transcript biotype "stringtie"
    7  remove ambiguous reference overlaps (7a, 7b, 7c)
    8  warnings and optional RefSeq flags
    9  merge with the reference and sort

Requires: polars, polars-bio.
"""
import argparse
import os
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference_gtf")
    parser.add_argument("assembly_gtf", help="gffcompare combined.gtf")
    parser.add_argument("tracking", help="gffcompare .tracking file")
    parser.add_argument("min_occurrence", type=float)
    parser.add_argument("min_tpm", type=float)
    parser.add_argument("output_prefix")
    parser.add_argument("platform", nargs="?", choices=["LR", "SR"], default="SR",
                        help="SR (default): long-read script; SR: short-read script (labels and single-exon rule)")
    parser.add_argument("threads", nargs="?", type=int, default=None,
                        help="threads for polars and polars-bio (default: all cores)")
    parser.add_argument("refseq_basename", nargs="?", default=None,
                        help="optional: <refseq_basename>.xr.gff and <refseq_basename>.nr.gff, only used for flags")
    return parser


def parse_arguments(argv=None):
    return build_parser().parse_args(argv)


# The arguments are parsed before polars is imported: polars fixes its thread pool size at import.
THREADS = parse_arguments(sys.argv[1:]).threads
if THREADS:
    os.environ["POLARS_MAX_THREADS"] = str(THREADS)

import polars as pl          # import after POLARS_MAX_THREADS is set
import polars_bio as pb      

if THREADS:
    # polars-bio runs its overlaps on DataFusion, which has its own number of parallel partitions
    pb.set_option("datafusion.execution.target_partitions", str(THREADS))
# GTF coordinates are 1-based and closed; pb.overlap() must treat every table that way.
pb.set_option(pb.POLARS_BIO_COORDINATE_SYSTEM_ZERO_BASED, False)

SCRIPT_VERSION, SCRIPT_DATE = "2.1", "26-01-2026"

# Chromosomes kept for the assembled transcripts: the autosomes plus X and Y.
MAIN_CHROMOSOMES = [str(number) for number in range(1, 23)] + ["X", "Y"]

# Attributes read from the reference GTF (gene_type only exists in GENCODE files).
REFERENCE_ATTRIBUTES = ["gene_id", "gene_name", "gene_biotype", "gene_type",
                        "transcript_id", "transcript_biotype", "exon_number"]
# Attributes gffcompare can write to combined.gtf.
GFFCOMPARE_ATTRIBUTES = ["transcript_id", "gene_id", "gene_name", "xloc", "ref_gene_id", "cmp_ref", "cmp_ref_gene",
                         "class_code", "tss_id", "num_samples", "contained_in", "oId", "exon_number"]

# Step 2: class codes whose structure is already in the reference (or that are artefacts).
ALREADY_ANNOTATED_CLASS_CODES = ["=", "c", "j", "m", "n", "e", "r", "s"]
# Step 3 correction: only these classes are real isoforms of the reference gene; x, i, y and p
# (antisense, intronic, containing, nearby) keep their own gene id.
SAME_GENE_CLASS_CODES = ["o", "k"]
# Step 3: attributes that gffcompare only writes on the transcript row.
METADATA_ON_TRANSCRIPT_ROW = ["gene_name", "oId", "cmp_ref", "class_code", "cmp_ref_gene", "ref_gene_id", "num_samples"]
# Step 7a: classes that should not overlap a reference exon.
NO_OVERLAP_CLASS_CODES = ["u", "x", "i", "y"]

# Columns of the output GTF (base_columns in R), in this order.
OUTPUT_COLUMNS = [
    "chrom", "source", "type", "start", "end", "score", "strand", "phase",
    "gene_id", "transcript_id", "exon_number", "gene_name", "gene_biotype", "transcript_biotype",
    "oId", "cmp_ref", "class_code", "tss_id", "num_samples", "num_samples_TPM_threshold",
    "cmp_ref_gene", "contained_in", "xloc", "transcript_evidence", "xr_overlap", "nr_overlap",
]

# The TPM is the third field from the end of a tracking entry: "...|<FPKM>|<TPM>|<coverage>|<length>".
TPM_FIELD_PATTERN = r"(?:^|\|)([^|]*)\|[^|]*\|[^|]*$"

# RefSeq accession -> Ensembl chromosome name (GRCh38).
REFSEQ_ACCESSION_TO_CHROMOSOME = {
    "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3", "NC_000004.12": "4",
    "NC_000005.10": "5", "NC_000006.12": "6", "NC_000007.14": "7", "NC_000008.11": "8",
    "NC_000009.12": "9", "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
    "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15", "NC_000016.10": "16",
    "NC_000017.11": "17", "NC_000018.10": "18", "NC_000019.10": "19", "NC_000020.11": "20",
    "NC_000021.9": "21", "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y",
}

# Sorting order of feature types within a transcript; unknown types come after these.
FEATURE_TYPE_PRIORITY = {"gene": 1, "transcript": 2, "exon": 3, "CDS": 4, "UTR": 4, "start_codon": 4}
STRANDED_FEATURE_TYPES = ["exon", "CDS", "UTR", "start_codon"]
FIXED_GTF_COLUMNS = ["chrom", "source", "type", "start", "end", "score", "strand", "phase"]


def log(message: str) -> None:
    """Progress messages go to stderr, so they end up in the Nextflow .command.err."""
    print(message, file=sys.stderr, flush=True)


# =====================================================================================================
# Reading input files 
# =====================================================================================================
def strand_or_star() -> pl.Expr:
    """'+' and '-' stay; '.' or a missing value becomes '*', Bioconductor's "no strand"."""
    return pl.when(pl.col("strand").is_in(["+", "-"])).then(pl.col("strand")).otherwise(pl.lit("*")).alias("strand")


def with_missing_columns(table: pl.DataFrame, names) -> pl.DataFrame:
    """Add every column in `names` that the table lacks, as an empty text column."""
    missing = [name for name in names if name not in table.columns]
    return table.with_columns([pl.lit(None, dtype=pl.String).alias(name) for name in missing])


def read_gtf_table(path, attributes: list) -> pl.DataFrame:
    """Read a GTF with polars-bio: 8 standard columns + one column per attribute.

    Attribute columns that are completely empty are dropped.
    """
    gtf = pb.scan_gtf(str(path), attr_fields=attributes)
    columns = gtf.collect_schema().names()
    missing = [name for name in FIXED_GTF_COLUMNS if name not in columns]
    if missing:
        raise KeyError(f"pb.scan_gtf() returned no column(s) {missing}; got {columns}")

    table = gtf.with_columns(pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64), strand_or_star()).collect()
    absent = [name for name in attributes
              if name in table.columns and table.get_column(name).null_count() == table.height]
    return table.drop(absent)


def use_ensembl_chromosome_names(table: pl.DataFrame) -> pl.DataFrame:
    """Rename UCSC-style chromosomes to Ensembl style: chr1 -> 1, chrX -> X, chrM -> MT."""
    without_chr_prefix = pl.col("chrom").str.replace(r"^chr", "")
    return table.with_columns(
        pl.when(without_chr_prefix == "M").then(pl.lit("MT")).otherwise(without_chr_prefix).alias("chrom")
    )


def sample_name_from_entry(entry: str) -> str:
    """'q1:sample1.03|sample1.03.1|…' -> 'sample1' (as extract_samples_from_tracking() in R)."""
    without_q_prefix = entry.split(":", 1)[-1]
    gene_in_sample = without_q_prefix.split("|", 1)[0]
    name, _, suffix = gene_in_sample.rpartition(".")
    return name if suffix.isdigit() else gene_in_sample


def read_tracking_table(path) -> pl.DataFrame:
    """A gffcompare .tracking file as one row per transcript with its TPM per sample.

    Columns: transcript_id, gene_id, ref_gene_id, ref_transcript_id, TPM_q1_<sample>, TPM_q2_<sample>, …
    """
    raw = pl.read_csv(path, separator="\t", has_header=False, quote_char=None, infer_schema=False)
    transcript_col, locus_col, reference_col, _class_col, *sample_cols = raw.columns

    reference_parts = pl.col(reference_col).str.splitn("|", 2)
    id_columns = [
        pl.col(transcript_col).alias("transcript_id"),
        pl.col(locus_col).alias("gene_id"),
        reference_parts.struct.field("field_0").alias("ref_gene_id"),
        reference_parts.struct.field("field_1").alias("ref_transcript_id"),
    ]
    tpm_columns = []
    for sample_number, column in enumerate(sample_cols, start=1):
        entries = raw.get_column(column)
        found = entries.filter(entries != "-")
        name = sample_name_from_entry(found[0]) if found.len() > 0 else None
        column_name = f"TPM_q{sample_number}_{name}" if name else f"TPM_q{sample_number}"
        tpm_columns.append(pl.col(column).str.extract(TPM_FIELD_PATTERN, 1).cast(pl.Float64).alias(column_name))
    return raw.select(id_columns + tpm_columns)


def read_refseq_exons(path) -> pl.DataFrame:
    """Exons from a RefSeq GFF3 file, with Ensembl-style chromosome names (main chromosomes only)."""
    features = (
        pb.scan_gff(str(path))
        .filter(pl.col("type") == "exon")
        .select(["chrom", "start", "end", "strand"])   # the polars-bio wrapper takes one list here
        .with_columns(pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64), strand_or_star())
        .collect()
    )
    chromosome_names = set(features["chrom"].unique().to_list())
    if chromosome_names & set(MAIN_CHROMOSOMES):
        pass
    elif any(name.startswith("chr") for name in chromosome_names):
        features = features.with_columns(pl.col("chrom").str.replace(r"^chr", ""))
    elif any(name.startswith("NC_") for name in chromosome_names):
        features = features.with_columns(pl.col("chrom").replace_strict(REFSEQ_ACCESSION_TO_CHROMOSOME, default=None))
    else:
        raise ValueError(f"Unknown chromosome names in {path}")
    return features.filter(pl.col("chrom").is_in(MAIN_CHROMOSOMES))


# =====================================================================================================
# Overlaps, sorting and writing
# =====================================================================================================
def overlapping_pairs(left: pl.DataFrame, right: pl.DataFrame, left_columns: list, right_columns: list,
                      same_strand: bool = True) -> pl.DataFrame:
    """All pairs of overlapping features, like GenomicRanges::findOverlaps(left, right).

    Returns the requested columns with the suffix "_left" / "_right", one row per overlapping pair.
    """
    result_columns = [f"{c}_left" for c in left_columns] + [f"{c}_right" for c in right_columns]
    if left.is_empty() or right.is_empty():
        return pl.DataFrame(schema={column: pl.String for column in result_columns})
    interval_columns = ["chrom", "start", "end", "strand"]
    pairs = pb.overlap(
        left.select(interval_columns + left_columns),
        right.select(interval_columns + right_columns),
        suffixes=("_left", "_right"),
        output_type="polars.DataFrame",
    )
    if same_strand:
        pairs = pairs.filter((pl.col("strand_left") == pl.col("strand_right"))
                             | (pl.col("strand_left") == "*") | (pl.col("strand_right") == "*"))
    return pairs.select(result_columns).unique()


def chromosome_sort_key() -> pl.Expr:
    """1..22, then X, Y, MT, then any other contig."""
    name = pl.col("chrom")
    return (
        pl.when(name.str.contains(r"^\d+$")).then(name.cast(pl.Int64, strict=False))
          .when(name == "X").then(pl.lit(23))
          .when(name == "Y").then(pl.lit(24))
          .when(name == "MT").then(pl.lit(25))
          .otherwise(pl.lit(1000)) # Anything else should have low priority
    )


def sort_gtf(table: pl.DataFrame) -> pl.DataFrame:
    """Sort a GTF table gene by gene."""
    # 1. genes: gene rows, plus the left-most transcript of genes without a gene row
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
    # 2. transcripts within each gene
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
    # 3. feature types; unknown types after the known ones, alphabetically
    alphabetical = sorted(table["type"].unique().to_list())
    rank_of_type = {feature: 4 + position for position, feature in enumerate(alphabetical, start=1)}
    rank_of_type.update(FEATURE_TYPE_PRIORITY)
    # 4. position within a transcript, 5' -> 3'
    position_in_transcript = (
        pl.when(pl.col("type").is_in(STRANDED_FEATURE_TYPES) & (pl.col("strand") == "-"))
          .then(-pl.col("start"))
          .otherwise(pl.col("start"))
    )
    # 5. sort; gene rows have no transcript rank (null) and come first, as with data.table
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

    gene_id and transcript_id come first; every other column with a value on that row follows.
    The lines are written straight to the file, without building the whole file as one string first.
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
        # one text column, written as is (no quoting, no header)
        lines.write_csv(handle, include_header=False, quote_style="never")


# =====================================================================================================
# The filter steps
# =====================================================================================================
def read_inputs(reference_path, assembly_path):
    """Step 0: the reference, and the assembled transcripts on the main chromosomes with a strand."""
    reference_gtf = use_ensembl_chromosome_names(read_gtf_table(reference_path, REFERENCE_ATTRIBUTES))
    # GENCODE: no gene_biotype, but gene_type -> use gene_type as gene_biotype
    if "gene_biotype" not in reference_gtf.columns and "gene_type" in reference_gtf.columns:
        reference_gtf = reference_gtf.rename({"gene_type": "gene_biotype"})
    reference_gtf = with_missing_columns(reference_gtf.drop("gene_type", strict=False), ["gene_biotype"])

    assembly_gtf = with_missing_columns(
        use_ensembl_chromosome_names(read_gtf_table(assembly_path, GFFCOMPARE_ATTRIBUTES)),
        ["cmp_ref", "class_code", "num_samples"],
    ).filter(pl.col("chrom").is_in(MAIN_CHROMOSOMES) & pl.col("strand").is_in(["+", "-"]))
    return reference_gtf, assembly_gtf


def label_reference_evidence(reference_gtf: pl.DataFrame, assembly_gtf: pl.DataFrame, prefix: str) -> pl.DataFrame:
    """Step 1: every reference transcript gets <prefix>_full_match / _partial_match / _other_match / _no_evidence."""
    best_evidence = (
        assembly_gtf
        .filter((pl.col("type") == "transcript") & pl.col("cmp_ref").is_not_null() & (pl.col("cmp_ref") != "."))
        .select("cmp_ref", "class_code")
        .unique()
        .group_by("cmp_ref")
        .agg(
            pl.when((pl.col("class_code") == "=").any()).then(pl.lit("full_match"))
              .when(pl.col("class_code").is_in(["j", "c"]).any()).then(pl.lit("partial_match"))
              .otherwise(pl.lit("other_match"))
              .alias("best_evidence")
        )
        .rename({"cmp_ref": "transcript_id"})
    )
    return (
        reference_gtf
        .join(best_evidence, on="transcript_id", how="left", maintain_order="left")
        .with_columns(
            pl.when(pl.col("type") != "transcript").then(pl.lit(None, dtype=pl.String))
              .when(pl.col("best_evidence").is_null()).then(pl.lit(f"{prefix}_no_evidence"))
              .otherwise(pl.lit(f"{prefix}_") + pl.col("best_evidence"))
              .alias("transcript_evidence")
        )
        .drop("best_evidence")
    )


def remove_already_annotated(assembly_gtf: pl.DataFrame):
    """Step 2: remove every row of transcripts whose class code says the structure is already known."""
    ids = (
        assembly_gtf
        .filter((pl.col("type") == "transcript") & pl.col("class_code").is_in(ALREADY_ANNOTATED_CLASS_CODES))
        .select("transcript_id")
        .unique()
    )
    return assembly_gtf.join(ids, on="transcript_id", how="anti"), ids


def attach_tracking(novel: pl.DataFrame, tracking_table: pl.DataFrame) -> pl.DataFrame:
    """Step 3: join the TPM per sample, copy transcript attributes to exon rows, and assign the gene id."""
    # every column name the two tables share is a join key
    join_keys = [column for column in tracking_table.columns if column in novel.columns]
    novel = novel.join(tracking_table, on=join_keys, how="left", maintain_order="left")
    novel = with_missing_columns(novel, ["gene_name", "ref_gene_id", "cmp_ref_gene"])

    novel = novel.with_columns([
        pl.col(column).fill_null(pl.col(column).filter(pl.col("type") == "transcript").first().over("transcript_id"))
        for column in METADATA_ON_TRANSCRIPT_ROW if column in novel.columns
    ])

    is_same_gene_class = pl.col("class_code").is_in(SAME_GENE_CLASS_CODES).fill_null(False)
    take_reference_gene = (pl.col("ref_gene_id") != "-") & pl.col("cmp_ref_gene").is_null() & is_same_gene_class
    novel = novel.with_columns(
        pl.when(take_reference_gene.is_null()).then(pl.lit(None, dtype=pl.String))
          .when(take_reference_gene).then(pl.col("ref_gene_id"))
          .otherwise(pl.col("gene_id"))
          .alias("gene_id"),
        # the reference gene's name would suggest that the other classes belong to it
        pl.when(is_same_gene_class).then(pl.col("gene_name")).otherwise(pl.lit(None, dtype=pl.String)).alias("gene_name"),
    )
    return novel.with_columns(pl.col("gene_name").fill_null(pl.col("gene_id")))


def remove_single_exon(novel: pl.DataFrame, only_without_reference_gene: bool):
    """Step 4: remove single-exon StringTie transcripts (LR: only those without a reference gene)."""
    exon_count = pl.len().over("transcript_id") - 1   # rows of the transcript minus its transcript row
    condition = (exon_count < 2) & (pl.col("source") == "StringTie")
    if only_without_reference_gene:
        condition = condition & (pl.col("ref_gene_id") == "-")
    ids = novel.filter(condition).select("transcript_id").unique()
    return novel.join(ids, on="transcript_id", how="anti"), ids


def filter_expression(novel: pl.DataFrame, min_occurrence: float, min_tpm: float):
    """Step 5: keep transcripts with TPM >= min_tpm in >= min_occurrence samples.

    Without a TPM threshold, num_samples >= min_occurrence (correction: R uses a strict >).
    """
    tpm_columns = [column for column in novel.columns if column.startswith("TPM_q")]
    if min_tpm > 0:
        samples_above_threshold = (
            pl.sum_horizontal([(pl.col(c).fill_null(0) >= min_tpm).cast(pl.Int32) for c in tpm_columns])
            if tpm_columns else pl.lit(0)
        )
        expressed = (
            novel
            .with_columns(samples_above_threshold.alias("num_samples_TPM_threshold"))
            .filter(pl.col("num_samples_TPM_threshold") >= min_occurrence)
        )
    elif min_occurrence > 0:
        expressed = novel.filter(pl.col("num_samples").cast(pl.Float64) >= min_occurrence)
    else:
        expressed = novel
    ids = novel.select("transcript_id").unique().join(expressed.select("transcript_id").unique(),
                                                      on="transcript_id", how="anti")
    return expressed, ids


def assign_biotype(novel: pl.DataFrame, reference_gtf: pl.DataFrame) -> pl.DataFrame:
    """Step 6: the biotype of the reference gene (first row of that gene), or 'stringtie' for new genes.

    Every added transcript gets transcript_biotype 'stringtie' (on all its rows): it is an assembled
    model, so its own biotype is unknown even when its gene is annotated.
    """
    biotype_per_gene = (
        reference_gtf
        .select("gene_id", pl.col("gene_biotype").alias("reference_biotype"))
        .unique(subset="gene_id", keep="first", maintain_order=True)
    )
    return (
        novel
        .drop("gene_biotype", "transcript_biotype", strict=False)
        .join(biotype_per_gene, on="gene_id", how="left", maintain_order="left")
        .with_columns(
            pl.col("reference_biotype").fill_null("stringtie").alias("gene_biotype"),
            pl.lit("stringtie").alias("transcript_biotype"),
        )
        .drop("reference_biotype")
    )


def overlap_filters(novel: pl.DataFrame, reference_gtf: pl.DataFrame):
    """Step 7: three strand-aware checks; a transcript flagged by any of them is removed.

    Returns the kept rows, and per check (7a, 7b, 7c, all) a table with the flagged transcript ids.
    """
    # 7a and 7b share one overlap: the reference genes each exon of a novel transcript overlaps
    exon_overlaps = overlapping_pairs(
        novel.filter(pl.col("type") == "exon"),
        reference_gtf.filter(pl.col("type") == "exon"),
        left_columns=["transcript_id", "class_code"], right_columns=["gene_id"],
    ).rename({"transcript_id_left": "transcript_id", "class_code_left": "class_code", "gene_id_right": "reference_gene"})

    # 7a: a "no overlap" class whose exon overlaps a reference exon anyway
    unexpected = (
        exon_overlaps.filter(pl.col("class_code").is_in(NO_OVERLAP_CLASS_CODES))
        .select("transcript_id").unique()
    )
    # 7b: exons overlapping exons of more than one reference gene
    multigene = (
        exon_overlaps.group_by("transcript_id")
        .agg(pl.col("reference_gene").n_unique().alias("genes_overlapped"))
        .filter(pl.col("genes_overlapped") > 1)
        .select("transcript_id")
    )
    # 7c: the whole span of a class i transcript against every reference feature, same strand
    same_strand_intronic = overlapping_pairs(
        novel.filter((pl.col("type") == "transcript") & (pl.col("class_code") == "i")), reference_gtf,
        left_columns=["transcript_id"], right_columns=["gene_id"],
    ).select(pl.col("transcript_id_left").alias("transcript_id")).unique()

    removed = pl.concat([unexpected, multigene, same_strand_intronic]).unique()
    kept = novel.join(removed, on="transcript_id", how="anti")
    return kept, unexpected, multigene, same_strand_intronic, removed


def refseq_flags(novel: pl.DataFrame, refseq_basename):
    """Step 8: flag transcripts with an exon overlapping a RefSeq XR / NR exon (nothing is removed)."""
    counts = {"xr": 0, "nr": 0}
    if not refseq_basename:
        return novel, counts
    novel_exons = novel.filter(pl.col("type") == "exon")
    for kind in ("xr", "nr"):
        hits = overlapping_pairs(novel_exons, read_refseq_exons(f"{refseq_basename}.{kind}.gff"),
                                 left_columns=["transcript_id"], right_columns=[])
        flagged = hits.select(
            pl.col("transcript_id_left").alias("transcript_id"), pl.lit(f"{kind}_hit").alias(f"{kind}_overlap"),
        ).unique()
        novel = (
            novel.join(flagged, on="transcript_id", how="left", maintain_order="left")
            .with_columns(pl.col(f"{kind}_overlap").fill_null("none"))
        )
        counts[kind] = novel.filter((pl.col("type") == "transcript") & (pl.col(f"{kind}_overlap") == f"{kind}_hit")).height
    return novel, counts


# =====================================================================================================
# Main
# =====================================================================================================
def main(argv=None) -> int:
    args = parse_arguments(argv)
    refseq_basename = (args.refseq_basename or "").strip() or None
    prefix = args.platform
    # long-read script: only single-exon transcripts without a reference gene; short-read script: all
    only_without_reference_gene = prefix == "LR"
    novel_label = f"{prefix}-novel"
    output_prefix = args.output_prefix

    log(f"polars uses {pl.thread_pool_size()} threads")
    log("step 0: reading the input files")
    reference_gtf, assembly_gtf = read_inputs(args.reference_gtf, args.assembly_gtf)
    log("step 1: reference evidence labels")
    reference_with_evidence = label_reference_evidence(reference_gtf, assembly_gtf, prefix)

    log("step 2: class code filter")
    novel, already_annotated_ids = remove_already_annotated(assembly_gtf)
    count_after_class_filter = novel["transcript_id"].n_unique()

    log("step 3: tracking file and gene ids")
    tracking_table = read_tracking_table(args.tracking)
    novel = attach_tracking(novel, tracking_table)

    log("step 4: single-exon filter")
    novel_multi_exonic, single_exon_ids = remove_single_exon(novel, only_without_reference_gene)

    log("step 5: expression filter")
    novel_expressed, low_expression_ids = filter_expression(novel_multi_exonic, args.min_occurrence, args.min_tpm)

    log("step 6: gene biotypes")
    novel_with_biotype = assign_biotype(novel_expressed, reference_gtf)

    log("step 7: overlap filters")
    novel_unambiguous, unexpected, multigene, intronic, _removed = overlap_filters(novel_with_biotype, reference_gtf)

    log("step 8: warnings and RefSeq flags")
    count_genes_on_both_strands = (
        novel_unambiguous.group_by("gene_id").agg(pl.col("strand").n_unique().alias("strands"))
        .filter(pl.col("strands") > 1).height
    )
    novel_flagged, refseq_counts = refseq_flags(novel_unambiguous, refseq_basename)

    # One entry per filter step: transcripts left after it, and its log lines. Drives the early exit and the log.
    filter_steps = [
        ("2 · class filter", count_after_class_filter,
         [("Novel transcripts - stranded, not in scaffolds:", count_after_class_filter)]),
        ("4 · single-exon filter", novel_multi_exonic["transcript_id"].n_unique(),
         [("Filtered monoexonic:", single_exon_ids.height)]),
        ("5 · expression filter", novel_expressed["transcript_id"].n_unique(),
         [("Filtered below expression and occurrence requirement:", low_expression_ids.height)]),
        ("7 · overlap filters", novel_unambiguous["transcript_id"].n_unique(),
         [("Filtered reference overlap (class u,x,i,y):", unexpected.height),
          ("Filtered reference multi-gene overlap:", multigene.height),
          ("Filtered class I same strand overlap:", intronic.height)]),
    ]
    early_exit_after = next((name for name, remaining, _ in filter_steps if remaining == 0), None)

    log("step 9: merge with the reference and sort")
    expression_table = novel_flagged.filter(pl.col("type") == "transcript").select(
        ["transcript_id", "gene_id"] + [c for c in novel_flagged.columns if c.startswith("TPM_q")])
    if early_exit_after is None:
        novel_labelled = novel_flagged.with_columns(pl.lit(novel_label).alias("transcript_evidence"))
        novel_transcripts_gtf = sort_gtf(novel_labelled.select([c for c in OUTPUT_COLUMNS if c in novel_labelled.columns]))
        extended_reference_gtf = (
            sort_gtf(pl.concat([novel_transcripts_gtf, reference_with_evidence], how="diagonal_relaxed"))
            .with_columns(pl.col("gene_name").fill_null(pl.col("gene_id")))
        )
        novel_ids = novel_transcripts_gtf["transcript_id"].unique().to_list()
    else:
        log(f"no transcripts remain after step {early_exit_after}: writing the reference only (as in R)")
        novel_transcripts_gtf = None
        extended_reference_gtf = reference_with_evidence
        novel_ids = []

    # ---- header and log --------------------------------------------------------------------
    gtf_header = [
        f"##GTF FILTERING {SCRIPT_VERSION} (last updated {SCRIPT_DATE} )",
        f"##reference_gtf={args.reference_gtf}",
        f"##assembled_transcriptome={args.assembly_gtf}",
        f"##tracking_file={args.tracking}",
        f"##refseq_gtf={refseq_basename or 'NA'}",
        f"##min_occurrence={args.min_occurrence:g}",
        f"##min_tpm={args.min_tpm:g}",
    ]

    log_lines = list(gtf_header)
    if early_exit_after == filter_steps[0][0]:
        log_lines.append("No novel transcripts remain after filtering.")
    else:
        for name, _, lines in filter_steps:
            log_lines += [f"{label} \t {value}" for label, value in lines]
            if name == early_exit_after:
                break
    log_lines.append("__________________________________________")
    if early_exit_after is None:
        log_lines += [f"TOTAL NOVEL TRANSCRIPTS ADDED TO REFERENCE: \t {len(novel_ids)}",
                      "__________________________________________"]
    evidence_counts = (
        extended_reference_gtf.filter(pl.col("type") == "transcript").drop_nulls("transcript_evidence")
        .group_by("transcript_evidence").len("transcripts").sort("transcript_evidence")
    )
    log_lines.append(f"Transcript evidence summary ({'long' if prefix == 'LR' else 'short'}-read support):")
    log_lines += [f"  - {row['transcript_evidence']:<18} : {row['transcripts']:>8}"
                  for row in evidence_counts.iter_rows(named=True)]
    if refseq_basename and early_exit_after is None:
        log_lines += [f"Info - Genes with transcripts on multiple strands: \t {count_genes_on_both_strands}",
                      f"Info - RefSeq xr overlap: \t {refseq_counts['xr']}",
                      f"Info - RefSeq nr overlap: \t {refseq_counts['nr']}"]

    # ---- outputs ----------------------------------------------------
    write_gtf(extended_reference_gtf, gtf_header, f"{output_prefix}.extended_reference.gtf")
    if novel_transcripts_gtf is not None:
        write_gtf(novel_transcripts_gtf, gtf_header, f"{output_prefix}.novel_transcripts.gtf")
    (expression_table if early_exit_after is None else expression_table.clear()).write_csv(
        f"{output_prefix}.tsv", separator="\t")
    Path(f"{output_prefix}.log").write_text("\n".join(log_lines) + "\n")

    log(f"done: {len(novel_ids)} novel transcripts added")
    return 0


if __name__ == "__main__":
    sys.exit(main())
