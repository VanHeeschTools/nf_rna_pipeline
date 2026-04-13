include { stringtie; stringtie_summary } from '../modules/stringtie'
include { mergeGTF; filterAnnotate; transcriptome_fasta } from '../modules/mergeTranscriptome'

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
    outdir             // Path to output directory

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
        stringtie(stringtie_input, chromosome_exclusion_list, reference_gtf, outdir)

        stringtie_summary(stringtie.out.stringtie_gtf.collect(),
                        reference_gtf,
                        outdir)

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
        mergeGTF(gtf_list, masked_fasta, reference_gtf, output_basename, outdir)

        gtf_novel = mergeGTF.out.merged_gtf
        gtf_tracking = mergeGTF.out.tracking

        refseq_input = refseq_gtf ? file("${refseq_gtf}*.gff") : []

        // Run filter annotate r script
        // TODO: Sort exons in gtf and add transcript biotype for stringtie tx
        filterAnnotate( reference_gtf,
                        refseq_input,
                        gtf_novel,
                        gtf_tracking,
                        min_occurrence,
                        min_tpm,
                        output_basename,
                        "${projectDir}/bin/",
                        outdir,
                        file("${projectDir}/bin/filter_annotate.R"),
                        file("${projectDir}/bin/filter_annotate_functions.R"))

        merged_filtered_gtf = filterAnnotate.out.gtf

        transcriptome_fasta(merged_filtered_gtf, masked_fasta, outdir)
        assembled_transcriptome_fasta = transcriptome_fasta.out
    } else {
        merged_filtered_gtf = null
        assembled_transcriptome_fasta = null
    }
    
    emit:
    stringtie_multiqc
    merged_filtered_gtf
    assembled_transcriptome_fasta
}
