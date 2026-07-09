include { indexLength; STAR } from '../modules/starAlign'

/*
* Runs star mapper with given input files, sort resulting bam files with samtools
*/
workflow ALIGN {
    take:
    reads              // Tuple, sample_id and location of input read files
    paired_end         // Bool, True if data is paired end
    reference_gtf      // Location of reference gtf
    star_index_basedir // Location of star index

    main:
    // Obtain read length to use correct star index
    index_length = indexLength(reads)

    // Map reads to genome using star
    STAR(reads,
        paired_end,
        index_length,
        reference_gtf,
        star_index_basedir
        )
        
    emit:
    bam             = STAR.out.sorted_bam           // Tuple of sample id and sorted bam file
    star_log        = STAR.out.star_log_final       // Log file of star for multiqc
    samtools_stats  = STAR.out.samtools_stats       // Statistics file of bam output for multiqc
}
