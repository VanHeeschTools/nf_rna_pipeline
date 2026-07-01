// Define process for stringtie
process stringtie {
    label "assembly"
    //publishDir "${outdir}/stringtie", mode: 'copy'

    input:
        tuple val(sample_id), val(stringtie), path(bam), val(paired) // Tuple of sample id and star result bam
        val chr                                         // Chromosome exclusion list 
        path reference_gtf                               // Reference gtf file location
        path outdir                                      // Path to output directory

    output:
        path "${sample_id}/${sample_id}.gff", emit: stringtie_gtf // Path to the output gtf of stringtie

    script:
        // Check if reads are unstranded and exit if they are
        if (stringtie == "unstranded") {
            println "Sample given to stringtie contains unstranded reads"
            // exit 1 
            strand = ""
        }else{
            strand = "--${stringtie}"
        }

        // Groovy function to add the chromosome exclusion list to the stringtie command
        def chromosome_exclusion = chr ? "-x ${chr}" : ""

        """
        stringtie ${bam[0]} \
            -G ${reference_gtf} \
            ${strand}\
            -M 0.45 \
            -a 9 \
            -m 50 \
            -f 0.05 \
            -g 40 \
            -j 2 \
            -s 99999 \
            -p ${task.cpus} \
            ${chromosome_exclusion} \
            -o "${sample_id}/${sample_id}.gff"
        """
}

process stringtie_summary {
    label "compareGTF"

    input:
        path gtf_list    
        path reference_gtf
        path outdir

    output:
        path "all_samples_stringtie_counts_mqc.tsv", emit: stringtie_multiqc

    script:
        """
        OUT="all_samples_stringtie_counts_mqc.tsv"
        echo -e "Sample\tTranscripts\tExons\tKnown_transcripts\tNovel_transcripts" >> "\$OUT"

        for GTF in ${gtf_list.join(' ')}; do
            # Extract sample_id from filename
            SAMPLE=\$(basename "\$GTF" .gtf)

            transcripts=\$(awk '\$3=="transcript"' "\$GTF" | wc -l)
            exons=\$(awk '\$3=="exon"' "\$GTF" | wc -l)

            gffcompare -r "$reference_gtf" -o "\${SAMPLE}_gffcmp" "\$GTF"
            ann="\${SAMPLE}_gffcmp.annotated.gtf"
            
            all=\$(grep -c \$'\ttranscript\t' "\$ann")
            known=\$(grep 'class_code "[=c]"' "\$ann" | grep -c \$'\ttranscript\t')
            novel=\$((all - known))

            echo -e "\$SAMPLE\t\$transcripts\t\$exons\t\$known\t\$novel" >> "\$OUT"
        done
        """
}
