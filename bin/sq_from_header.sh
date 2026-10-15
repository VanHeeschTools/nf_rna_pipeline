#!/usr/bin/env bash
# Reads SAM header lines on stdin and prints one line per @SQ record with
# the sequence name, M5 checksum (lowercase) and reference URI, tab separated.
# Used to match CRAM files to the reference fasta they were compressed against.

awk -F'\t' '
$1 == "@SQ" {
    sn = ""; m5 = ""; ur = ""
    for (i = 2; i <= NF; i++) {
        tag = substr($i, 1, 2)
        value = substr($i, 4)
        if (tag == "SN") sn = value
        else if (tag == "M5") m5 = tolower(value)
        else if (tag == "UR") ur = value
    }
    print sn "\t" m5 "\t" ur
}'
