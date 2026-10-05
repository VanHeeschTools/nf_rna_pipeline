// Define process to list the sequence names and M5 checksums of a CRAM header
process cramHeaderSq {
    label "samtools"
    tag "${meta.id}"

    input:
        tuple val(meta), path(cram) // Tuple of sample meta map and cram file

    output:
        tuple val(meta), path("${meta.id}.sq.tsv"), emit: sq // Table of sequence name, M5 checksum and reference URI

    script:
        """
        set -o pipefail
        samtools view -H "${cram}" | sq_from_header.sh > "${meta.id}.sq.tsv"
        """

    stub:
        // Stub runs use plain text test files containing only @SQ header lines
        """
        grep '^@SQ' "${cram}" | sq_from_header.sh > "${meta.id}.sq.tsv"
        """
}

// Define process to list the sequence names and M5 checksums of a reference fasta
process referenceSq {
    label "samtools"
    tag "${fasta_file.name}"

    input:
        tuple val(fasta), path(fasta_file), path(fai) // Tuple of original fasta path, fasta file and its index

    output:
        tuple val(fasta), path("reference.sq.tsv"), emit: sq // Table of sequence name, M5 checksum and reference URI

    script:
        """
        set -o pipefail
        samtools dict "${fasta_file}" | sq_from_header.sh > reference.sq.tsv
        """

    stub:
        // Stub runs use plain text test files containing only @SQ header lines
        """
        grep '^@SQ' "${fasta_file}" | sq_from_header.sh > reference.sq.tsv
        """
}

// Define process to convert a CRAM file to paired fastq files
process cramToFastq {
    label "samtools"
    tag "${meta.id}"

    input:
        tuple val(meta), path(cram), path(fasta), path(fai) // Tuple of sample meta map, cram file and matching reference

    output:
        tuple val(meta.id), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: fastq // Tuple of sample id and paired fastq files
        path "${meta.id}_{singleton,other}.fastq.gz", emit: unpaired // Reads without a mate, not used downstream

    script:
        def prefix = meta.id
        """
        set -o pipefail

        # Group mates together first, CRAM files are usually sorted by coordinate
        # and samtools fastq only keeps R1 and R2 in sync for name-grouped input
        samtools collate \
            -u -O \
            -@ ${task.cpus} \
            --reference "${fasta}" \
            "${cram}" \
            collate_tmp \
        | samtools fastq \
            -@ ${task.cpus} \
            -1 "${prefix}_R1.fastq.gz" \
            -2 "${prefix}_R2.fastq.gz" \
            -s "${prefix}_singleton.fastq.gz" \
            -0 "${prefix}_other.fastq.gz" \
            -

        # Only paired-end CRAM input is supported
        if [ -z "\$(gzip -cd "${prefix}_R2.fastq.gz" | head -c 1)" ]; then
            echo "ERROR: no paired reads found in ${cram}, single-end CRAM input is not supported" >&2
            exit 1
        fi
        """

    stub:
        def prefix = meta.id
        """
        for read_file in R1 R2 singleton other; do
            printf '@${prefix}_read\\nACGT\\n+\\nIIII\\n' | gzip > "${prefix}_\${read_file}.fastq.gz"
        done
        """
}
