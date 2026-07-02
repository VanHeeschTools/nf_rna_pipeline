process starAlignChimeric {
    label "fusions"

    input:
        tuple val(sample_id), path(reads)
        val paired_end
        path reference
        val outdir

    output:
        tuple val("${sample_id}"), path("${sample_id}.*.bam"), emit: bam
        path "${sample_id}.*"

    script:
        if (paired_end == true){
            """
            STAR --runThreadN $task.cpus \
            --genomeDir ${reference}/STAR_index* \
            --genomeLoad NoSharedMemory \
            --readFilesIn "${reads[0]}" "${reads[1]}" \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outFileNamePrefix "${sample_id}." \
            --outBAMcompression 0 \
            --outFilterMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimSegmentMin 10 \
            --chimOutType WithinBAM HardClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --outTmpDir /tmp/star_tmp/

            """
        } else {
            error  "`--fusions` not supported when `--paired_end false`"
        }
}

process runArriba{
    label "fusions"

    input:
        tuple val(sample_id), path(bam), path(vcf)
        path fa
        path fai
        path gtf
        path blacklist
        path whitelist
        path protein_domains
        path outdir

    output:
        path "${sample_id}_fusions*", emit: fusions

    script:
        // If blacklist is [], this evaluates to false
        def blacklist_options = blacklist ? "-b ${blacklist}" : "-f blacklist"
        def whitelist_options = whitelist ? "-k ${whitelist} -t ${whitelist}" : ""
        def domains_options   = protein_domains ? "-p ${protein_domains}" : ""
        def wgs_options       = vcf ? "-d ${vcf}" : ""

        """
        mkdir -p ${sample_id}/
        /arriba_v2.4.0/arriba -x ${bam} \
        -o ${sample_id}_fusions.tsv \
        -O ${sample_id}_fusions.discarded.tsv \
        -a ${fa} \
        -g ${gtf} \
        ${wgs_options} ${blacklist_options} ${whitelist_options} ${domains_options}
        """
}


