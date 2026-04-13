include { validateParameters; paramsSummaryLog; samplesheetToList; } from 'plugin/nf-schema'
include { checkInputFiles; is_paired_end } from "../modules/helperFunctions.nf"
include { QC } from '../subworkflows/QC'
include { ALIGN } from '../subworkflows/ALIGN'
include { ASSEMBLY } from '../subworkflows/ASSEMBLY'
include { FUSIONS } from '../subworkflows/FUSIONS'
include { EXPRESSION } from '../subworkflows/EXPRESSION'
include { MULTIQC } from '../modules/multiqc.nf'

workflow RNASEQ {

    // Initialise workflow
    // Set module toggles based on run_mode
    run_qc         = params.run_mode == "qc" ? true  : (params.run_mode == "qc_restart" ? false : params.qc)
    run_align      = params.run_mode == "qc" ? false : params.align 
    run_assembly   = params.run_mode == "qc" ? false : params.assembly
    run_merge      = params.run_mode == "qc" ? false : params.merge
    run_fusions    = params.run_mode == "qc" ? false : params.fusions
    run_expression = params.run_mode == "qc" ? false : params.expression

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Validate files
    checkInputFiles()

    // Check if input is paired_end or single end and prevent mixing
    boolean paired_end_check = is_paired_end(params.input)
    // Turn into channel to give to modules
    paired_end = channel.from(paired_end_check).first()

    // Create input channel from samplesheet
    ch_input = channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    // Create channel with sample id and RNA reads (single or paired)
    reads = ch_input
        .filter { meta, _filename_1, _filename_2 ->
            meta.filetype == "fastq" && meta.sequence == "rna"
        }
        .map { meta, fastq_1, fastq_2 ->
            def read_files = fastq_2 ? [fastq_1, fastq_2] : [fastq_1]
            tuple(meta.id, read_files)
        }

    // Build patient_id maps for wgs matching
    rna_map = ch_input
        .filter { meta, _filename_1, _filename_2 -> meta.sequence == "rna" }
        .map { meta, _filename_1, _filename_2 -> [meta.subject, meta.id] }

    multiqc_files = channel.empty()

    //TODO handle bamfiles from samplesheet

    //Handle WGS vcfs from samplesheet
    vcfs = ch_input
        .filter { meta, _filename_1, _filename_2 ->
            meta.filetype == "vcf" && meta.sequence == "dna"
        }
        .map { meta, filename_1, _filename_2 ->
            tuple(meta.subject, filename_1)
        }

    // Match to rna sample id
    // Note: this will only work if there's one vcf per patient
    // If multiple vcfs per patient, vcf can me merged beforhand or match to the right rna sample by adjusting patient ids (i.e. patient1_diagnostic, patient1_relapse for both rna and wgs) 

    vcfs = vcfs
        .join(rna_map, remainder: true)    // If no vcf returns [sample, null]
        .map { _patient, vcf, rna_sample_id ->
            tuple(rna_sample_id, vcf)
        }

    /*
    * Step 01: QC
    */
    if (run_qc) {
        QC(reads,
        paired_end_check, 
        paired_end,
        params.kallisto_index,
        params.reference_gtf,
        params.strandedness_check,
        params.outdir,
        params.store_trimmed_reads)
        
        star_input = QC.out.trimmed_reads
        strand = QC.out.strandedness
        
        multiqc_files = multiqc_files.mix(QC.out.fastp_json)
        multiqc_files = multiqc_files.mix(QC.out.strand_file)

    } else {
        // Look for trimmed reads from previous QC run
        trimmed_glob = "${params.outdir}/**/*{R1,R2}*_trimmed.{fastq.gz,fq.gz}"
        
        // Checl if file(trimmed_glob) expected to be list of files is not empty
        if (!file(trimmed_glob).isEmpty()) {
            if (!params.bam_files) {
                log.info "Using trimmed reads from previous QC: ${trimmed_glob}"
                star_input = channel.fromFilePairs(trimmed_glob, size: paired_end ? 2 : 1, checkIfExists: true)
            } else {
                log.info "BAM files provided, skipping use of trimmed reads."
                star_input = null
            }
        } else {
            log.info "No trimmed reads found. Using raw input reads."
            star_input = reads
        }

        // Look for strandedness summary file
        default_strand_file = "${params.outdir}/check_strandedness/strandedness_all.txt"
        if (params.strand_info) {

            log.info "Using strandedness info from user-provided file: ${params.strand_info}"

            channel.fromPath(params.strand_info, checkIfExists: true)
            .splitText { line -> line.trim() }

        } else if (params.run_mode == "qc_restart" && !file(default_strand_file).isEmpty()) {

            log.info "Using strandedness from previous QC: ${default_strand_file}"

            channel.fromPath(default_strand_file, checkIfExists: true)
                .splitCsv(sep: "\t", header: false)
                .set { strand }
        } else {

            log.info "No strandedness provided or found, setting strand to `null`. WARNING: It will conflict with assembly step."
            log.info "You can specify path to strandedness info file using --strand_info parameter."

            strand = null
        }
    }

    /*
    *Step 02: Align
    */
    if (run_align){
        ALIGN(star_input, paired_end, params.reference_gtf, params.star_index_basedir, params.outdir)
        bam = ALIGN.out.bam
        multiqc_files = multiqc_files.mix(ALIGN.out.star_log)
        multiqc_files = multiqc_files.mix(ALIGN.out.samtools_stats)

    } else {
        // Look for alignment files in case star has been run previously
        default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
        if (params.bam_files) {
            bam = channel.fromFilePairs(params.bam_files, size: 1, checkIfExists: true)
        } else if (!file(default_bams).isEmpty()) {
            bam = channel.fromFilePairs(default_bams, size: 1, checkIfExists: true)
        } else {
            bam = null
        }
    }

    /*
    * Step 03: Assemble transcriptome
    */
    if ((run_assembly || run_merge) && paired_end_check == true) {

        // Join the strand info with the bam file to prevent sample mixing
        stringtie_input = strand.join(bam)
        .map { row ->
                def sample_id = row[0]    // Sample Id
                def strand_label = row[1] // Strand value
                def file_path = row[3]    // BAM file path
                [sample_id, strand_label, file_path, paired_end_check]
            }

        ASSEMBLY(stringtie_input,
                params.sample_gtf_list,
                params.reference_gtf,
                params.refseq_gtf,
                params.chr_exclusion_list,
                params.masked_fasta,
                params.output_basename,
                params.min_occurrence,
                params.min_tpm,
                params.outdir)

        assembled_gtf = ASSEMBLY.out.merged_filtered_gtf
        assembled_fasta = ASSEMBLY.out.assembled_transcriptome_fasta
        multiqc_files = multiqc_files.mix(ASSEMBLY.out.stringtie_multiqc)
    } else {
        assembled_gtf = null
        assembled_fasta = null
    }
    
    /*
    * Step 04: Fusion calling
    */

    // Merge star input with wgs vcfs to prevent sample mixing
    if (run_fusions && paired_end_check == true) {
        FUSIONS(star_input,
            vcfs,
            paired_end,
            params.arriba_reference,
            params.outdir)
    }

    /*
    * Step 05: Expression
    */
    if (run_expression) {
        // Join the strand info with the bam file to prevent sample mixing
        if (run_align){
            featurecounts_input = strand.join(bam)
            .map { row ->
                    def sample_id = row[0]    // Sample Id
                    def strand_code = row[2]  // Strand code value
                    def file_path = row[3]    // BAM file path
                    [sample_id, strand_code, paired_end_check, file_path]
                }
        } else {
            featurecounts_input = null
            println "Skipping featurecounts as no new alignments were made."
        }
        
        EXPRESSION(star_input,
            featurecounts_input,
            assembled_gtf,
            assembled_fasta,
            params.expression_mode,
            paired_end_check,
            params.reference_transcriptome,
            params.reference_gtf,
            params.output_basename,
            params.outdir)
        multiqc_files = multiqc_files.mix(EXPRESSION.out.salmon_multiqc)
    }

    /*
    * Step 06: Immune landscape
    */
    
    /*
    * Step 07: MultiQC
    */
    multiqc_config_file = file("${projectDir}/${params.multiqc_config}")
    MULTIQC (multiqc_files.collect(), multiqc_config_file, params.outdir) 
}

