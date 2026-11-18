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
 * Checks the samplesheet for ids that consist only of digits (e.g. "123" or "0123").
 * The samplesheet validation reads them as numbers instead of text, which fails with an
 * unclear message and would drop leading zeros. Stops with a clear message before anything runs.
 *
 * @param samplesheet path to the CSV samplesheet
 */
def check_numeric_ids(samplesheet) {
    def lines = file(samplesheet).readLines().findAll { line -> line.trim() }
    if (!lines) {
        return
    }

    def header = lines[0].split(',', -1)*.trim()
    def id_columns = ["sample_id", "subject_id", "group_id"].findAll { column -> column in header }
    def numeric_ids = []

    lines.drop(1).eachWithIndex { line, index ->
        def values = line.split(',', -1)*.trim()
        id_columns.each { column ->
            def column_index = header.indexOf(column)
            def value = column_index < values.size() ? values[column_index] : ""
            if (value ==~ /-?\d+(\.\d+)?/) {
                numeric_ids << [line: index + 2, column: column, value: value]
            }
        }
    }

    if (numeric_ids) {
        def message = ["", "The samplesheet contains ids that consist only of digits, which the pipeline can't use:"]
        message += numeric_ids.collect { entry -> "  line ${entry.line}: ${entry.column} '${entry.value}'".toString() }
        message << "Add at least one letter to these ids (for example 'S${numeric_ids[0].value}') and use the same ids in related files such as strand_info."
        message << "No analysis was started."
        error message.join('\n')
    }
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
    // Check strandedness: QC determines it, otherwise it comes from strand_info or a previous QC run.
    // Assembly and merge need it (together with bams) for every sample.
    def default_strand =  "${params.outdir}/check_strandedness/strandedness_all.txt"
    def run_qc = params.run_mode == "qc_restart" ? false : (params.run_mode == "qc" || params.qc)
    def needs_strand = params.run_mode != "qc" && !run_qc && (params.align || params.assembly || params.merge)
    if (needs_strand) {
        if (params.strand_info) {
            if (!file(params.strand_info, type : "file").exists()) {
                error "--strand_info: file doesn't exist, check path ${params.strand_info}\nNo analysis was started."
            }
            log.info "strand info   : ${params.strand_info}"
        } else if (file(default_strand, type : "file").exists()) {
            log.info "strand info   : ${default_strand}"
        } else {
            def steps = ["align", "assembly", "merge"].findAll { step -> params[step] }.join(', ')
            error """
                No strandedness information available, which is needed for: ${steps}.
                QC is not run (qc = false or run_mode = qc_restart) and ${default_strand} doesn't exist.
                Either run QC (qc = true), or give a strandedness file with `strand_info`
                (one line per sample: 'sample_id rf', 'sample_id fr' or 'sample_id unstranded').
                No analysis was started.
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
    // Assembly and merge need bams together with the strandedness of each sample
    if ((params.assembly || params.merge) && params.run_mode != "qc" && !bam_avail) {
        error  """
            Assembly and merge need bam files, but none were found in ${default_bams}
            If `align = false`, define your bam location with `bam_files=[path_to_file]` in `params.config`
            No analysis was started.
            """.stripIndent()
    }

    if (params.assembly){

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

    // Check the long-read transcriptome
    if (params.longread_gtf) {
        check_files("longread_gtf", params.longread_gtf, "file")
        if (!params.merge_longread) {
            log.info "longread gtf  : given, but merge_longread = false, the transcriptome is built from short reads only"
        } else if (!params.merge) {
            log.info "longread gtf  : given, but merge = false, so it is not merged with a short-read transcriptome"
        } else {
            log.info "longread gtf  : ${params.longread_gtf} (merged with the short-read transcriptome, end_tolerance ${params.end_tolerance} nt)"
        }
    }

    // Check HLA typing parameters
    if (params.hla_typing && params.run_mode != "qc") {
        if (!params.align && !bam_avail) {
            error  """
                HLA typing needs bam files, but none were found in ${default_bams}
                If `align = false`, define your bam location with `bam_files=[path_to_file]` in `params.config`
                """.stripIndent()
        }
        // arcasHLA names its outputs after the part of the file name before the first dot
        def ids_with_dot = samples
            .findAll { row -> row[0].sequence == "rna" && row[0].id.toString().contains('.') }
            .collect { row -> row[0].id }
        if (ids_with_dot) {
            error "HLA typing does not support sample ids containing a dot ('.'), rename these samples: ${ids_with_dot.join(', ')}"
        }
        log.info "hla typing    : arcasHLA genotype and expression"
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

/*
 * Reads one line of a strandedness file (strand_info or check_strandedness/strandedness_all.txt):
 * "sample_id<whitespace>strand_type" with strand_type rf, fr or unstranded.
 * Returns the same tuple as getStrandtype: [sample_id, strand_type, featurecounts code].
 */
def parse_strand_line(String line) {
    def strand_codes = ["rf": 1, "fr": 2, "unstranded": 0]
    def fields = line.trim().split(/\s+/)
    def strand_type = fields.size() > 1 ? fields[1].toLowerCase() : ""
    if (fields.size() < 2 || !strand_codes.containsKey(strand_type)) {
        error """
            Invalid line in the strandedness file: '${line.trim()}'
            Expected 'sample_id strand_type' with strand_type 'rf', 'fr' or 'unstranded'.
            """.stripIndent()
    }
    return tuple(fields[0], strand_type, strand_codes[strand_type])
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
