include { CRAM_TO_FASTQ } from './CRAM_TO_FASTQ'

/*
* Collects the RNA-seq reads from the samplesheet.
* Fastq rows are used as is, CRAM rows are converted to paired fastq files first.
*/
workflow PREPARE_READS {
    take:
    ch_input              // Tuple of sample meta map, filename_1 and filename_2 from the samplesheet
    conversion_references // Path to text file listing candidate reference fastas for CRAM input, or null

    main:
    rna_input = ch_input
        .filter { meta, _filename_1, _filename_2 -> meta.sequence == "rna" }

    // Fastq input: tuple of sample id and single or paired read files
    fastq_reads = rna_input
        .filter { meta, _filename_1, _filename_2 -> meta.filetype == "fastq" }
        .map { meta, fastq_1, fastq_2 ->
            def read_files = fastq_2 ? [fastq_1, fastq_2] : [fastq_1]
            tuple(meta.id, read_files)
        }

    // CRAM input: converted to paired fastq files
    crams = rna_input
        .filter { meta, _filename_1, _filename_2 -> meta.filetype == "cram" }
        .map { meta, cram, _filename_2 -> tuple(meta, cram) }

    CRAM_TO_FASTQ(crams, conversion_references)

    emit:
    reads = fastq_reads.mix(CRAM_TO_FASTQ.out.reads) // Tuple of sample id and input read file(s)
}
