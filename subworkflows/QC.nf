include { fastp }         from '../modules/fastp'
include { checkStrand }   from '../modules/strandedness'
include { getStrandtype } from '../modules/helperFunctions.nf'

//Run fastp on the input reads and obatains their strandedness
workflow QC {
    take:
    reads                     // Input fastq files
    paired_end_check          // Bool, true if paired end data
    paired_end                // Channel bool, true if paired end data
    kallisto_index            // Path to the kallisto index file
    reference_gtf             // Path to the reference gtf 
    strandedness_check        // Val, number of reads to use for strandedness check
    store_trimmed_reads       // Bool, true if trimmed reads need to be stored

    main:
    // Run fastp and sets the output fastq file to variable trimmed_reads
    fastp(reads, paired_end, store_trimmed_reads)

    // Run strandedness
    if (paired_end_check == true){
        strand = checkStrand(reads, kallisto_index, reference_gtf, strandedness_check).strand
        
        // Map strand info to strandedness format used downstream
        strandedness = strand.map{ it -> getStrandtype(it) }
        strand_file = checkStrand.out.strand_file

        // Collect strandedness info into a single file
        strandedness
        .map { tup -> "${tup[0]}\t${tup[1]}" }
        .collectFile(
            name: 'strandedness_all.txt',
            storeDir: "${params.outdir}/check_strandedness/",
            newLine: true,
            sort: true
        )
    } else {
        strand = null
        strandedness = null
        strand_file = null        
    }

    emit:
    strandedness  // File, strandedness of the data
    trimmed_reads   = fastp.out.fastq_files   // Tuple, sample and trimmed fastq files
    fastp_json      = fastp.out.fastp_json    // Fastp stats
    strand_file   // Text file containing the strandeness of all input samples
}