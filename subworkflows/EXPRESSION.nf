include { salmon_index; salmon_quasi; salmon_tables; featurecounts} from '../modules/salmon'

/*
* Run salmon steps to obtain expression statistics
*/
workflow EXPRESSION {

    take:
    reads           // Tuple of sample id and input read file(s)
    featurecounts_input
    assembled_gtf   // Path to assembled transcriptome gtf file
    assembled_fasta // Path to assembled transcriptome sequences file
    mode            // Expression mode to run
    paired_end      // Bool, is data paired end or not
    transcriptome   // Path to the input transcriptome file
    reference_gtf   // Path to the input reference gtf file
    output_basename // File prefix given to the salmon_tables results
    outdir          // Path to output directory

    main:

    // Initialize outputs as empty
    salmon_multiqc = channel.empty()
    
    // Use stringtie created transcriptome if set to true, otherwise use reference
    if (params.created_transcriptome_expression && assembled_gtf != null){
        input_gtf = assembled_gtf.first()
        input_transcriptome = assembled_fasta.first()
        println "Using assembled transcriptome for quantification"
    } else {
        input_gtf = reference_gtf
        input_transcriptome = transcriptome
    }
    
    if (mode =~ /.*sq.*/) {
        // Create the salmon index of the given transcriptome
        salmon_index(input_transcriptome)

        // Run salmon and write paths of quant.sf output files to a text file
        salmon_quasi(reads, paired_end, salmon_index.out, outdir)


        // Write the paths of the salmon_quasi output files to a text file
        ch_quants = salmon_quasi.out.quant.collect()

        // Run the salmon_tables Rscript to obtain expression tables
        salmon_tables(ch_quants, input_gtf, output_basename, outdir)
        salmon_multiqc = salmon_tables.out.salmon_multiqc
    }


    // Run featurecounts if bam files for input exist and strand info is available
    
    if (mode =~ /.*fc.*/){
        print("Running featurecounts for expression quantification")
        featurecounts(featurecounts_input, input_gtf, outdir)
    }

    emit:
    salmon_multiqc
}
