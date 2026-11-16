include { stringtie; stringtie_summary } from '../modules/stringtie'
include { mergeGTF; filterAnnotate; transcriptome_fasta; select_lr_novel_records; trmap; lr_sr_combine_transcriptome; } from '../modules/mergeTranscriptome'

workflow ASSEMBLY {
    take:
    stringtie_input
    sample_gtf_list    // List of previously created gtf files to be merged
    reference_gtf      // Path to input reference gtf file
    refseq_gtf         // Path to input refseq gtf file
    chr_exclusion_list // Path to chromosome exclustion list
    masked_fasta       // Path to input reference fasta file
    output_basename    // Val containing the id/name given to the output files
    min_occurrence     // Val contatining the minimum occurence of transcripts for filtering
    min_tpm            // Val containing the minium tpm of transcripts for filtering
    longread_gtf       // Path to custom gtf output from the nf_longread pipeline
    end_tolerance      // Val max nt difference at which LR and SR transcript ends still count as the same

    main:
    // Run stringtie unless paths to precomputed individual sample GTF are provided
    if (!sample_gtf_list & params.assembly) {

        // Locate chromosome exclusion list
        if(chr_exclusion_list) {
        // Groovy command to join the chr exclusion list on ","
            chromosome_exclusion_list = file("${projectDir}/${chr_exclusion_list}").readLines().join(",")
        } else {
            chromosome_exclusion_list = null
        }

        // Run stringtie
        stringtie(stringtie_input, chromosome_exclusion_list, reference_gtf)

        stringtie_summary(stringtie.out.stringtie_gtf.collect(),
                        reference_gtf)

        stringtie_multiqc = stringtie_summary.out.stringtie_multiqc

        // Collect GTF files and create a list file
        gtf_list = stringtie.out.stringtie_gtf.collect()
    } else {
        stringtie_multiqc = null
    }

    // Merges the gtf files created by stringtie
    if (params.merge) {

        // Load gtf list file if not null
        if (sample_gtf_list) {
            gtf_list = channel.fromPath("${sample_gtf_list}")
            gtf_list
            .splitText{ line -> line.trim() }
            .take(1)
            .ifEmpty { error "Could not find sample GTF files in: ${sample_gtf_list}" }
        }

        // Run merge process
        mergeGTF(gtf_list, masked_fasta, reference_gtf, output_basename)

        gtf_merged = mergeGTF.out.merged_gtf
        gtf_tracking = mergeGTF.out.tracking

        refseq_input = refseq_gtf ? file("${refseq_gtf}*.gff") : []

        // Run filter annotate Python script
        if (longread_gtf) {
            filter_annotate_min_occurence = 1
        } else{
            filter_annotate_min_occurence = min_occurrence
        }

        filterAnnotate( reference_gtf,
                        refseq_input,
                        gtf_merged,
                        gtf_tracking,
                        filter_annotate_min_occurence,
                        min_tpm,
                        output_basename)

        merged_filtered_gtf = filterAnnotate.out.gtf

        if (longread_gtf){
            // Filter long read gtf to obtain only novel transcripts
            select_lr_novel_records(longread_gtf)
            // Overlap long read novel transcripts to short read novel transcripts
            trmap(filterAnnotate.out.gtf_novel, select_lr_novel_records.out.lr_novel_gtf)
            // Combine long-read and short-read transcriptoem
            lr_sr_combine_transcriptome(
                filterAnnotate.out.gtf,
                params.longread_gtf,
                trmap.out.trmap_tab_file,
                min_occurrence,
                end_tolerance
            )
            transcriptome_fasta(lr_sr_combine_transcriptome.out.gtf, masked_fasta)
            output_gtf = lr_sr_combine_transcriptome.out.gtf
        } else{
            transcriptome_fasta(merged_filtered_gtf, masked_fasta)
            output_gtf = merged_filtered_gtf
        }

        assembled_transcriptome_fasta = transcriptome_fasta.out
    } else {
        merged_filtered_gtf = null
        assembled_transcriptome_fasta = null
    }
    
    emit:
    stringtie_multiqc
    output_gtf
    assembled_transcriptome_fasta
}
