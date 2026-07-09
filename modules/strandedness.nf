// Define process for checking strandedness
process checkStrand {
    label "qc"

    input:
        tuple val(sample_id), path(reads) // Tuple containing sample id and read paths
        path kallisto_index                // Path to kallisto index dir
        path reference_gtf                 // Path to input reference gtf file
        val strandedness_check     // Val, number of reads to use for strandedness check

    output:
        tuple val("${sample_id}"), stdout, emit: strand
        path "${sample_id}.txt", emit: strand_file

    script:
        """
        check_strandedness \
            -g ${reference_gtf} \
            -n ${strandedness_check} \
            -r1 ${reads[0]} \
            -r2 ${reads[1]} \
            -k "${kallisto_index}" >> "${sample_id}.txt"
        strandedness=\$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print \$NF}')
        printf "%s\t%s" "${sample_id}" "\$strandedness"
        rm -r stranded_test_* >&2
        """
}
