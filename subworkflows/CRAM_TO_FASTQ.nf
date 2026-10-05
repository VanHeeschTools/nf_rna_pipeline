include { cramHeaderSq; referenceSq; cramToFastq } from '../modules/cramToFastq'
include { select_cram_reference } from '../modules/helperFunctions.nf'

/*
* Converts RNA-seq CRAM files to paired fastq files.
* The reference fasta of each CRAM comes from the samplesheet `reference` column, or is
* selected from `cram_references` by matching sequence names and M5 checksums of the CRAM header.
*/
workflow CRAM_TO_FASTQ {
    take:
    crams            // Tuple of sample meta map and cram file
    cram_references  // Path to text file listing candidate reference fastas (one per line), or null

    main:
    // Collect candidate reference fastas from cram_references and the samplesheet
    if (cram_references) {
        listed_references = channel.fromPath(cram_references, checkIfExists: true)
            .splitText()
            .map { line -> line.trim() }
            .filter { line -> line && !line.startsWith('#') }
    } else {
        listed_references = channel.empty()
    }

    samplesheet_references = crams
        .filter { meta, _cram -> meta.reference }
        .map { meta, _cram -> meta.reference.toString() }

    // Use full paths so the same fasta given in both places is only checked once
    reference_fastas = listed_references
        .mix(samplesheet_references)
        .map { fasta -> file(fasta).toUriString() }
        .unique()
        .map { fasta -> tuple(fasta, file(fasta), file("${fasta}.fai")) }

    // Obtain sequence names and M5 checksums of the CRAM headers and the reference fastas
    cramHeaderSq(crams)
    referenceSq(reference_fastas)

    // Gather the tables of all reference fastas into a single item: [[fasta, lines], ...]
    reference_sq = referenceSq.out.sq
        .map { fasta, sq_table -> [fasta, sq_table.readLines()] }
        .toList()
        .map { references -> [references] }

    // Select the reference fasta matching each CRAM
    cram_with_reference = crams
        .map { meta, cram -> tuple(meta.id, meta, cram) }
        .join(cramHeaderSq.out.sq.map { meta, sq_table -> tuple(meta.id, sq_table.readLines()) })
        .combine(reference_sq)
        .map { _sample_id, meta, cram, cram_sq, references ->
            def fasta = select_cram_reference(meta, cram_sq, references)
            log.info "cram input    : ${meta.id} uses reference ${fasta}"
            tuple(meta, cram, file(fasta), file("${fasta}.fai"))
        }

    // Convert CRAM to paired fastq files
    cramToFastq(cram_with_reference)

    emit:
    reads = cramToFastq.out.fastq.map { sample_id, fastq_1, fastq_2 -> tuple(sample_id, [fastq_1, fastq_2]) } // Tuple of sample id and [R1, R2] fastq files
}
