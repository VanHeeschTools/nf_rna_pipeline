#!/usr/bin/env python3
"""
Selects the HLA genes that arcasHLA customize and quant can quantify for one sample.

arcasHLA customize builds the sample specific HLA reference from the cDNA sequences in its fixed
dat/ref/cDNA.json (default grouping 'protein-group'). Typed alleles that are missing there get no
sequence: when both alleles of a gene are missing, arcasHLA quant fails with a KeyError, and when
one allele is missing it is silently left out of the quantification. Only genes whose typed
alleles are all present are quantified; the other genes are written to a table.

Usage: arcashla_quant_genes.py <sample.genotype.json> <arcasHLA installation folder> <skipped genes tsv>
Prints the genes to quantify as a comma separated list (empty when none).
"""
import json
import os
import sys

genotype_json, arcas_dir, skipped_tsv = sys.argv[1:4]

# Use the 2-field conversion of arcasHLA itself, as done by arcasHLA customize
sys.path.insert(0, os.path.join(arcas_dir, "scripts"))
from arcas_utilities import process_allele  # noqa: E402

with open(os.path.join(arcas_dir, "dat", "ref", "cDNA.json")) as handle:
    cdna = json.load(handle)
with open(genotype_json) as handle:
    genotype = json.load(handle)

genes = []
skipped = []
for gene, alleles in sorted(genotype.items()):
    missing = sorted({process_allele(allele, 2) for allele in alleles if not cdna.get(process_allele(allele, 2))})
    if missing:
        skipped.append((gene, ",".join(alleles), ",".join(missing)))
    else:
        genes.append(gene)

with open(skipped_tsv, "w") as handle:
    handle.write("gene\ttyped_alleles\tmissing_in_arcashla_cdna\treason\n")
    for gene, alleles, missing in skipped:
        handle.write(f"{gene}\t{alleles}\t{missing}\tno cDNA sequence in arcasHLA dat/ref/cDNA.json, gene not quantified\n")

print(",".join(genes))
