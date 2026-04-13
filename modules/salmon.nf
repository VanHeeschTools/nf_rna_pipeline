// Create Salmon index file
process salmon_index {
    label "salmon"

    input:
        path transcriptome // Existing transcriptome file

    output:
        path 'salmon_index'

    script:
        """
        salmon index --threads $task.cpus -t $transcriptome -i salmon_index -k 31
        """
}

// Run Salmon quant mode
process salmon_quasi {
    label "salmon"
    publishDir "${outdir}/salmon", mode: 'copy', pattern: "${sample_id}/quant.sf"

    input:
        tuple val(sample_id), path(reads)  // Tuple of sample id and input read file(s)
        val paired_end                     // Bool, is data paired end or not
        path salmon_index                  // Path to the salmon index
        val outdir                         // Path to the output directory

    output:
        //path "${sample_id}/quant.sf", emit: quant
        path("${sample_id}_quant.sf"), emit: quant

    script:
        if (paired_end == true){
            // Set salmon quant paired end input arguments
            quant_input = """-1 "${reads[0]}" -2 "${reads[1]}" """
        } else {
            // Set salmon quant single end input arguments
            quant_input = """-r "${reads[0]}" """
        }
        """
        salmon quant \
        --libType "A" \
        --validateMappings \
        --gcBias \
        --quiet \
        --numGibbsSamples 30 \
        --threads $task.cpus \
        -i "${salmon_index}" \
        ${quant_input} \
        --output "${sample_id}"

        # Rename the specific file for uniqueness
        mv ${sample_id}/quant.sf ./${sample_id}_quant.sf
        """
}

// Create satistics table using the salmon output
process salmon_tables {
    label "salmon_tables"
    publishDir "${outdir}/salmon", mode: 'copy', pattern: "${prefix}*"

    input:
        path quants, stageAs: "quants/*"
        path gtf
        val prefix
        val outdir

    output:
        path "${prefix}*"
        path "${prefix}_multiqc_summary_mqc.tsv", emit: salmon_multiqc

    script:
        def path_list = quants.collect { "\$(pwd)/${it}" }.join('\n') 
        """
        echo "${path_list}" > quant_paths.txt
 
        salmon_cohort_tables.R \
        quant_paths.txt \
        ${gtf} \
        ${prefix}
        """
}

// Run featurecounts 
process featurecounts {
    label "featurecounts"
    publishDir "${outdir}/featurecounts", mode: 'copy'

    input:
        tuple val(sample_id), val(strand), val(paired_end), path(bam) // Tuple of sample id and input read file(s)
        path reference_gtf                                   // Path to the input reference gtf file
        path outdir                                          // Path to output directory

    output:
        path "${sample_id}/*"

    /*
    optional parameters to be considered:
    -t : specify feature type in GTF annotation (exon by default)
    -g : specify attribute type in GTF annotation (gene_id by default)
    */

    script:
        paired_flag = (paired_end == true) ? "-p" : ""

        """
        mkdir -p ${sample_id}
        featureCounts \
        ${paired_flag} \
        -t gene \
        -s ${strand} \
        -a ${reference_gtf} \
        -o ${sample_id}/featurecounts_result.out \
        -T $task.cpus \
        ${bam}
        """
}
