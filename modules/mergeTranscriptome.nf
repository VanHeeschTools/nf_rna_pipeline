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

        # If the combined file wasn't created, find the single evaluated GTF and rename it
        if [ ! -f "${output_basename}.combined.gtf" ]; then
            echo "Single file detected. Creating fallback combined.gtf..."
            cp ${output_basename}.*.gtf ${output_basename}.combined.gtf
        fi
        """
}


// Define process for transcript filtering and annotation
process filterAnnotate {
    label "assembly"

    input:
        path reference_gtf  // Path to the input reference gtf file
        path refseq_files   // Path to input refseq gtf file
        path gtf_novel      // Path to the merged gtf file
        path gtf_tracking   // Path to the tracking file created by the merge step
        val min_occurrence  // Val contatining the minimum occurence of transcripts for filtering
        val min_tpm         // Val containing the minium tpm of transcripts for filtering
        val output_basename // Val containing output basename
        path scripts_dir    // Path location of input R scripts
        path "filter_annotate.R"
        path "filter_annotate_functions.R"

    output:
        path "${output_basename}.extended_reference.gtf", emit: gtf
        path "${output_basename}.novel_transcripts.gtf", emit: gtf_novel
        path "${output_basename}.log"
        path "${output_basename}.tsv"

    when:
        task.ext.when == null || task.ext.when

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
        "${output_basename}" \
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


process select_lr_novel_records {
    label "salmon_tables"

    input:
        path longread_gtf

    output:
        path "lr_novel.transcript_records.gtf", emit: lr_novel_gtf

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    library(rtracklayer)

    gtf <- import("${longread_gtf}")
    keep_ids <- unique(gtf\$transcript_id[!is.na(gtf\$transcript_evidence) &
                                        gtf\$transcript_evidence == "LR_novel"])
    export(gtf[gtf\$transcript_id %in% keep_ids], "lr_novel.transcript_records.gtf")
    """
}

process trmap{
    label "compareGTF"

    input:
        path merged_filtered_novel_gtf
        path lr_novel_gtf
    
    output:
        path "trmap_sr_vs_lr.out", emit: trmap_out_file
        path "trmap_sr_vs_lr.tab", emit: trmap_tab_file

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        trmap -o "trmap_sr_vs_lr.out" \
                "${merged_filtered_novel_gtf}" \
                "${lr_novel_gtf}"

        trmap -T -o "trmap_sr_vs_lr.tab" \
                "${merged_filtered_novel_gtf}" \
                "${lr_novel_gtf}"
        """
}

process lr_sr_combine_transcriptome{
    debug true  
    label "salmon_tables"

    input:
        path sr_gtf_file
        path lr_gtf_file
        path trmap_out_file
        path trmap_tab_file
        val min_occurrence
        path "filter_annotate_functions.R"
    
    output:
        path "unified_lr_sr_novel.gtf"
        path "unified_lr_sr.gtf"
    
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        lr_sr_combine_transcriptomes_trmap.R \
            "${sr_gtf_file}" \
            "${lr_gtf_file}" \
            "${trmap_out_file}" \
            "${trmap_tab_file}" \
            "." \
            "${min_occurrence}" \
        """

}