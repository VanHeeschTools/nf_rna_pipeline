include { arcasReference; arcasExtract; arcasGenotype; arcasQuant; arcasMerge } from '../modules/arcashla'

/*
* Immune landscape of the samples.
* HLA typing and HLA allele expression with arcasHLA from the STAR alignments.
* TODO: CIBERSORTx and IPASS
*/
workflow IMMUNE {
    take:
    bam             // Tuple of sample id and sorted bam file
    paired_end      // Bool, true if data is paired end
    output_basename // Val, prefix of the cohort HLA tables

    main:
    // Check the IMGT/HLA database first, so a missing database stops the run before the long extraction
    arcasReference()

    // Add the bam index when it is next to the bam, arcasHLA indexes the bam itself otherwise
    hla_input = bam
        .map { sample_id, bam_files ->
            def bam_file = [bam_files].flatten()[0]
            def bai_file = file("${bam_file}.bai")
            tuple(sample_id, bam_file, bai_file.exists() ? bai_file : [])
        }
        .combine(arcasReference.out.versions)
        .map { sample_id, bam_file, bai_file, _versions -> tuple(sample_id, bam_file, bai_file) }

    // Extract chromosome 6 and unmapped reads and predict the HLA genotype
    arcasExtract(hla_input, paired_end)
    arcasGenotype(arcasExtract.out.fastq, paired_end)

    // Combine the genotypes of all samples into cohort tables
    merge_input = arcasGenotype.out.genotype
        .map { _sample_id, genotype -> genotype }
        .mix(arcasGenotype.out.genes)
        .collect()
    arcasMerge(merge_input, "${output_basename}.hla")

    // Quantify expression of the typed HLA alleles
    // Samples without a predicted genotype (e.g. too few HLA reads) are skipped, arcasHLA customize fails on them
    quant_input = arcasExtract.out.fastq
        .join(arcasGenotype.out.genotype)
        .filter { sample_id, _fastq, genotype ->
            def typed = !(genotype.text.trim() in ["", "{}"])
            if (!typed) {
                log.warn "HLA typing: no HLA genotype predicted for ${sample_id}, HLA expression is skipped for this sample"
            }
            typed
        }
    arcasQuant(quant_input, paired_end)

    emit:
    hla_genotypes = arcasMerge.out.genotypes // Table with the HLA genotype of all samples
}
