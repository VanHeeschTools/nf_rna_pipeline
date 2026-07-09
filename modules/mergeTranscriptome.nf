// Define process for GTF merging
process mergeGTF {
    label "compareGTF"

    input:
        path gtf_list         // List of previously created gtf files to be merged
        path masked_fasta      // Path to input reference fasta file
        path reference_gtf     // Path to input reference gtf file
        val output_basename   // Val containing the id/name given to the output files

    output:
        path "${output_basename}*"
        path "${output_basename}.combined.gtf", emit: merged_gtf
        path "${output_basename}.tracking", emit: tracking

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        ls *.gff > gtflist.txt
        
        gffcompare \
            -V \
            -r ${reference_gtf} \
            -s ${masked_fasta} \
            -o "${output_basename}" \
            -i gtflist.txt

        """
}


// Define process for transcript filtering and annotation
process filterAnnotate {
    label "assembly"

    input:
        path reference_gtf   // Path to the input reference gtf file
        path refseq_files      // Path to input refseq gtf file
        path gtf_novel      // Path to the merged gtf file
        path gtf_tracking   // Path to the tracking file created by the merge step
        val min_occurrence  // Val contatining the minimum occurence of transcripts for filtering
        val min_tpm         // Val containing the minium tpm of transcripts for filtering
        val output_basename // Val containing output basename
        path scripts_dir     // Path location of input R scripts
        path "filter_annotate.R"
        path "filter_annotate_functions.R"

    output:
        path "${output_basename}_novel_filtered.gtf", emit: gtf
        path "${output_basename}_novel_filtered.log"
        path "${output_basename}_novel_filtered.tsv"

    script:
	def refseq_prefix = refseq_files ? refseq_files[0].name.replace(".xr.gff", "").replace(".nr.gff", "") : ""
        def refseq_arg = refseq_prefix ? "\"${refseq_prefix}\"" : ""
        """
        filter_annotate.R \
        "${reference_gtf}" \
        "${gtf_novel}" \
        "${gtf_tracking}" \
        "${min_occurrence}" \
        "${min_tpm}" \
        "${output_basename}_novel_filtered" \
        "${scripts_dir}" \
        ${refseq_arg}
        """
}

// Creates a fasta file of the transcript sequence using the reference fasta file and the transcriptome gtf
process transcriptome_fasta {
    label "gffread"

    input:
        path merged_filtered_gtf // Merged and filtered transcriptome file
        path masked_fasta        // Path to input reference fasta file

    output:
        file "stringtie_transcriptome.fa"

    script:
        """
        gffread -w stringtie_transcriptome.fa -g ${masked_fasta} ${merged_filtered_gtf}
        """
}
