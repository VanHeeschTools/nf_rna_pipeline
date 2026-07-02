// Define process to assign best index length
process indexLength{
    label "alignment"

    input:
        tuple val(sample_id), path(reads) // Tuple containing sample id and read paths
    
    output:
        stdout
    
    script:
        """
        # Check first 10k reads for read length for star index
        read_length=\$(gunzip -c "${reads[0]}" | \
        head -n 10000 | \
        awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count - 1}')

        # Make sure read length picks are correct
        if [[ "\${read_length}" =~ "." ]]; then
        read_length=\${read_length%.*}
        fi

        if (( \${read_length} >= 70 && \${read_length} <= 85 )); then
        used_index=74
        elif (( \${read_length} >= 86 && \${read_length} <= 111 )); then
        used_index=99
        elif (( \${read_length} >= 112 && \${read_length} <= 137 )); then
        used_index=124
        elif (( \${read_length} >= 138 && \${read_length} <= 163 )); then
        used_index=149
        else
        used_index=99
        fi

        printf "%s" "\${used_index}"
        """

}

// Define process for alignment with STAR
process STAR {
    label "alignment"

    input: 
        tuple val(sample_id), path(reads) // Tuple containing sample id and read paths
        val paired_end                    // Bool set to true if input reads are paired end
        val usedIndex                     // Val showing read length for the star index
        path reference_gtf                 // Input reference gtf file
        path star_index_basedir            // Path to star index dir
        path outdir                        // Path to output directory

    output:
        //path "${sample_id}/${sample_id}.*"
        //tuple val("${sample_id}"), path("${sample_id}/${sample_id}.*.bam"), emit: bam
        tuple val(sample_id), path("${sample_id}/${sample_id}*.Aligned.sortedByCoord.out.bam"), emit:sorted_bam
        path "${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam.bai", emit:bam_bai
        path "${sample_id}/${sample_id}.Log.final.out", emit:star_log_final
        path "${sample_id}/${sample_id}_stats.txt", emit:samtools_stats
    
    script:
        // STAR params defined as concatenated strings singe triple quote multiline declaration generates newline
        def star_params = "--readFilesCommand zcat " +
                        "--twopassMode Basic --runDirPerm All_RWX " +
                        "--outFilterType BySJout --outSAMunmapped Within " +
                        "--outSAMattributes NH HI AS nM NM MD jM jI MC ch " +
                        "--outSAMstrandField intronMotif --outSAMtype BAM Unsorted " +
                        "--outFilterMismatchNmax 6 --outTmpKeep None " +
                        "--alignSJoverhangMin 10 --outFilterMultimapNmax 10 " +
                        "--outFilterScoreMinOverLread 0.75" 

        def bam = "${sample_id}.Aligned.out.bam"
        def new_bam = "${sample_id}.Aligned.sortedByCoord.out.bam"

        if (paired_end == true){
            // Set paired-end input parameter
            star_input = """--readFilesIn "${reads[0]}" "${reads[1]}" """
        } else {
            // Set single-end input parameter
            star_input = """--readFilesIn "${reads}" """
        }
        
        """
        mkdir -p ${sample_id}
        mkdir -p tmp/

        # Use STAR for mapping the reads
        time STAR --genomeDir "${star_index_basedir}/${usedIndex}nt" \
        --sjdbGTFfile ${reference_gtf} \
        ${star_input} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id} \
        --outFileNamePrefix "${sample_id}/${sample_id}." \
        --runThreadN $task.cpus ${star_params} \
        --outTmpDir /tmp/star_tmp/


        # Sort BAM
        time samtools sort \
        -@ $task.cpus  \
        -l 9 \
        -o "${sample_id}/${new_bam}" \
        -T "tmp/" \
        "${sample_id}/${bam}"

        # Delete tmp dir and unsorted bam file
        rm -r tmp/
        rm ${sample_id}/${bam} 

        # Create mapping statistics with samtools
        samtools stats -@ $task.cpus "${sample_id}/${new_bam}" > "${sample_id}/${sample_id}_stats.txt"

        # Index the bam with samtools
        samtools index -@ $task.cpus "${sample_id}/${new_bam}"
        """
}
