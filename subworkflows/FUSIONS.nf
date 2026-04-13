include { starAlignChimeric; runArriba} from '../modules/arriba'

workflow FUSIONS {
    take:
    reads             // Tuple: (sample_id, [fastq1, fastq2])
    vcf               // Tuple: (sample_id, vcf) or (sample_id, null) if no vcf available
    paired_end        // Channel bool, paired end or not
    arriba_reference  // Path arriba reference location
    outdir            // Path to output dir

    main:
    // Location to the Arriba reference
    def ref_files = file("${arriba_reference}/**")

    // Check if both fasta and index exists for Arriba reference 
    def fa  = ref_files.find { it.name.endsWith('viral.fa') }
    def fai  = ref_files.find { it.name.endsWith('.fai') }
    //fa_with_index = [ fa, file("${fa}.fai") ]
    def gtf = ref_files.find { it.name.endsWith('.gtf') }

    // Take optional files and passing empty file if empty
    def blacklist       = ref_files.find { it.name.contains('blacklist') } ?: []
    def whitelist       = ref_files.find { it.name.contains('known_fusions') } ?: []
    def protein_domains = ref_files.find { it.name.contains('protein_domains') } ?: []

    ch_vcf_cleaned = vcf.map { it -> 
        def meta = it[0]
        def vcf_file = it[1] ?: [] 
        return [meta, vcf_file]
    }

    // Run star mapper
    starAlignChimeric(reads, paired_end, arriba_reference, outdir)
    arriba_input = starAlignChimeric.out.bam.join(ch_vcf_cleaned)

    // Run Arriba
    runArriba(arriba_input, 
        fa,
        fai, 
        gtf,
        blacklist, 
        whitelist, 
        protein_domains, 
        outdir
    )
}
