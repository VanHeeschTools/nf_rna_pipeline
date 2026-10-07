def printHeader () {
    def logMessage =  """
        .-------------------------------------------------------.
        |                _____                 _      _     _   |
        | _ _ ___ ___   |  |  |___ ___ ___ ___| |_   | |___| |_ |
        || | | .'|   |  |     | -_| -_|_ -|  _|   |  | | .'| . ||
        | \\_/|__,|_|_|  |__|__|___|___|___|___|_|_|  |_|__,|___||
        '-------------------------------------------------------'

        ${workflow.manifest.name} ${workflow.manifest.version}
        ==========================
        """
        log.info logMessage.stripIndent()
}


/*
 * Finds sample ids that occur more than once in the samplesheet.
 *
 * @param samples list of [meta, filename_1, filename_2] rows from samplesheetToList
 * @return list of duplicated sample ids, empty if all are unique
 */
def duplicate_sample_ids(List samples) {
    return samples
        .collect { row -> row[0].id }
        .countBy { sample_id -> sample_id }
        .findAll { _sample_id, count -> count > 1 }
        .keySet()
        .toList()
}

/*
 * Lists the reference fastas given for CRAM input: the samplesheet `reference`
 * column and the fastas listed in `params.conversion_references`.
 *
 * @param samples list of [meta, filename_1, filename_2] rows from samplesheetToList
 * @return unique list of reference fasta paths as Strings
 */
def conversion_reference_fastas(List samples) {
    def fastas = samples
        .findAll { row -> row[0].filetype == "cram" && row[0].conversion_reference }
        .collect { row -> row[0].conversion_reference.toString() }

    if (params.conversion_references) {
        fastas += file(params.conversion_references).readLines()
            .collect { line -> line.trim() }
            .findAll { line -> line && !line.startsWith('#') }
    }

    return fastas.unique()
}

/*
 * Turns "name<TAB>m5<TAB>uri" lines into a set of "name<TAB>m5" keys.
 */
def sq_keys(List lines) {
    return lines.collect { line ->
        def fields = line.split('\t', -1)
        "${fields[0]}\t${fields[1]}".toString()
    } as Set
}

/*
 * Selects the reference fasta of a CRAM file by comparing the sequence names and
 * M5 checksums of the CRAM header with those of the candidate reference fastas.
 * When the samplesheet gives a reference, only that fasta is checked.
 *
 * @param meta       sample meta map
 * @param cram_sq    lines "name<TAB>m5<TAB>uri" for every @SQ record of the CRAM header
 * @param references list of [fasta, lines] with the same table for each candidate fasta
 * @return path of the matching reference fasta as String
 */
def select_conversion_reference(meta, List cram_sq, List references) {
    def cram_fields = cram_sq.collect { line -> line.split('\t', -1) }
    def cram_uris = cram_fields.collect { fields -> fields.size() > 2 ? fields[2] : "" }.findAll { uri -> uri }.unique()
    def given_reference = meta.conversion_reference ? file(meta.conversion_reference.toString()).toUriString() : null

    // Without M5 checksums the reference can't be verified
    if (!cram_fields || cram_fields.any { fields -> fields.size() < 2 || !fields[1] }) {
        if (given_reference) {
            log.warn "CRAM header of ${meta.id} has no M5 checksums, using ${given_reference} from the samplesheet without verification"
            return given_reference
        }
        error """
            CRAM header of ${meta.id} has no M5 checksums, so its reference can't be matched automatically.
            Fill the `conversion_reference` column in the samplesheet for this sample.
            Reference used to create the CRAM: ${cram_uris ? cram_uris.join(', ') : 'unknown'}
            """.stripIndent()
    }

    def cram_keys = sq_keys(cram_sq)
    def candidates = given_reference ? references.findAll { fasta, _lines -> fasta == given_reference } : references
    def match = candidates.find { _fasta, lines -> sq_keys(lines).containsAll(cram_keys) }

    if (!match) {
        def reason = given_reference
            ? "The reference given in the samplesheet (${given_reference}) does not match the CRAM header."
            : "None of the fastas in `conversion_references` matches the CRAM header."
        error """
            No matching reference fasta found for CRAM sample ${meta.id}.
            ${reason}
            Sequence names and M5 checksums of the reference must match all @SQ lines of the CRAM header.
            Reference used to create the CRAM: ${cram_uris ? cram_uris.join(', ') : 'unknown'}
            """.stripIndent()
    }

    return match[0]
}

// Checks existence of given file path, can be a file or a directory
def check_files(name, path, type) {
    if(!path) {
        error "When running merge without assembly you must provide `${name}`."
    } else {
        def file_to_check = null
        if (type == "dir") {
            file_to_check = file(path, type: "dir")
        } else {
            file_to_check = file(path, type: "file")
        }
        if (file_to_check != null) {
            // Check if it's a list of files
            if (file_to_check instanceof List) {
                // Check each file in the list
                file_to_check.flatten().each { file ->
                    if (!file.exists()) {
                        error "--${name}: ${type} doesn't exist, check path ${file}"
                    }
                }
            } else {
                // Check the single file
                if (!file_to_check.exists()) {
                    error "--${name}: ${type} doesn't exist, check path ${path}"
                }
            }
        } else {
            error "--${name}: No files found at ${path}"
        }
    }
}

// Checks existence of files listed in this function
// samples: list of [meta, filename_1, filename_2] rows from samplesheetToList
def checkInputFiles(List samples) {
    // Check inputs
    def default_strand =  "${params.outdir}/check_strandedness/strandedness_all.txt"
    if (!params.qc && ( params.align || params.assembly )) {
        if (!params.strand_info && file(default_strand, type : "file").exists()) {
            log.info "strand info   : ${default_strand}".stripIndent()
        } else if (params.strand_info && file(params.strand_info, type : "file").exists()) {
            log.info "strand info   : ${params.strand_info}".stripIndent()
        } else {
            error  """
                No strand information found in ${default_strand}
                If `qc = false`, define your strandedness info file with `strand_info` in `params.config`
                """.stripIndent()
        }
    }

    // Locate bams
    def default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
    def bam_avail = true
    if (!params.align) {
        if (!params.bam_files && !file(default_bams).isEmpty()) {
            log.info "bam files     : ${default_bams}".stripIndent()
        } else if (params.bam_files && !file(params.bam_files).isEmpty()) {
            log.info "bam files     : ${params.bam_files}".stripIndent()
        } else {
            bam_avail = null
        }
    }

    // Check reference_gtf
    check_files("reference_gtf", params.reference_gtf, "file")

    // Check kallisto index
    if (params.qc){
        check_files("kallisto_index", params.kallisto_index, "file")
    }

    // Check star_index_basedir
    if (params.align){
        check_files("star_index_basedir", params.star_index_basedir, "dir")
    }

    // Check bams and references for assembly steps
    if (params.assembly){
        if (!bam_avail) {
            error  """
                No bam found in ${default_bams}
                If `align = false`, define your bam location with `bam_files=[path_to_file]` in `params.config`
                """.stripIndent()
        }

        check_files("refseq_gtf", "${params.refseq_gtf}*", "file")
        check_files("masked_fasta", "${params.masked_fasta}", "file")
    }

    // Check sample gtf list for merge without assembly
    if (!params.assembly && params.merge) {
        if(!params.sample_gtf_list) {
            error "When running merge without assembly you must provide `sample_gtf_list`.".stripIndent()
        } else {
            def sample_gtf_list = file(params.sample_gtf_list, type: "file")
            if (!sample_gtf_list.exists()) {
                error  "--sample_gtf_list: file doesn't exist, check path ${params.sample_gtf_list}".stripIndent()
            }
        }
    }

    // Check references for build_annotaiton
    if (params.build_annotation) {
        check_files("twobit", "${params.twobit}*", "file")
    }

    // Check fusion parameters
    if (params.fusions){
        check_files("arriba_reference", params.arriba_reference, "dir")
        check_files("arriba_reference STAR_index_*", "${params.arriba_reference}STAR_index*", "dir")
        check_files("arriba_reference *.gtf", "${params.arriba_reference}*.gtf", "file")
        check_files("arriba_reference *.fa", "${params.arriba_reference}*.fa", "file")
    }

    // Check expression parameters
    if (params.expression) {
        assert params.expression_mode in ["sq", "fc", "sqfc"]: "`expression_mode` must be one of the following: sq, fc, sqfc"

        if (!bam_avail & (params.expression_mode != "sq")) {
            log.info "expression mode  : ${params.expression_mode} --> sq [forced to quasi-mapping, no bam avail]"
        }
    }

    // Check conversion reference fastas for CRAM input (only used to convert CRAM to fastq)
    def rna_crams = samples.findAll { row -> row[0].sequence == "rna" && row[0].filetype == "cram" }
    if (rna_crams) {
        def without_reference = rna_crams.findAll { row -> !row[0].conversion_reference }.collect { row -> row[0].id }
        if (without_reference && !params.conversion_references) {
            error """
                No conversion reference fasta available for CRAM samples: ${without_reference.join(', ')}
                Fill the `conversion_reference` column in the samplesheet, or set `conversion_references` to a text file listing candidate reference fastas (one path per line).
                """.stripIndent()
        }
        if (params.conversion_references) {
            check_files("conversion_references", params.conversion_references, "file")
        }
        conversion_reference_fastas(samples).each { fasta ->
            check_files("conversion reference fasta", fasta, "file")
            if (!file("${fasta}.fai").exists()) {
                error "--conversion reference fasta: index ${fasta}.fai doesn't exist, create it with `samtools faidx ${fasta}`"
            }
        }
        // The CRAM index (.crai) is not needed, filename_2 is ignored for CRAM rows
        def with_filename_2 = rna_crams.findAll { row -> row[2] }.collect { row -> row[0].id }
        if (with_filename_2) {
            log.warn "filename_2 is ignored for CRAM input (no .crai index needed), samples: ${with_filename_2.join(', ')}"
        }
        log.info "cram input    : ${rna_crams.size()} sample(s), converted to fastq before QC"
    }

    // Check for duplicate sample ids
    def duplicates = duplicate_sample_ids(samples)
    if (duplicates) {
        error "\nERROR: Duplicate sample_ids found in the input samplesheet: ${duplicates.join(', ')}\n"
    }

    log.info "\n==========================\n"
}

/*
 * Checks whether all RNA samples in a samplesheet are paired-end or single-end.
 * CRAM samples are converted to paired-end fastq files, so they count as paired-end.
 * Throws an error if there's a mix.
 *
 * @param samples list of [meta, filename_1, filename_2] rows from samplesheetToList
 * @return true if all samples are paired-end, false if all are single-end
 */
def is_paired_end(List samples) {
    def paired_flags = samples
        .findAll { row -> row[0].sequence == "rna" && row[0].filetype in ["fastq", "cram"] }
        .collect { row -> row[0].filetype == "cram" ? true : (row[2] ? true : false) }
        .unique()

    if (paired_flags.isEmpty()) {
        error "Samplesheet contains no RNA fastq or cram rows."
    }
    if (paired_flags.size() > 1) {
        error "Mixed paired-end and single-end data detected in samplesheet! Please only use one type of data. Note: CRAM input counts as paired-end."
    }

    return paired_flags[0]
}

// Groovy function to check if strand type exists and to set the strand type if found
def getStrandtype(strand_tuple) {
    def (sample_id, strand) = strand_tuple
    def strand_type = null

    if (strand ==~ /.*RF\/fr-firststrand/) {
        strand_type = tuple(sample_id, "rf", 1)
    } else if (strand ==~ /.*FR\/fr-secondstrand/) {
        strand_type = tuple(sample_id, "fr", 2)
    } else {
        strand_type = tuple(sample_id, "unstranded", 0)
        // Optionally warn or exit here
        // println "WARNING: Data for ${sample_id} is unstranded"
        // exit 1
    }

    return strand_type
}
