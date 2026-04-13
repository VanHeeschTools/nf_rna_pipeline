process MULTIQC {
    label "multiqc"
    publishDir "${outdir}/multiqc/", mode: 'copy'

    input:
        path multiqc_files, stageAs: "?/*"
        path multiqc_config
        path outdir

    output:
        path "*multiqc_report.html", emit: multiqc_report

    script:
        """
        multiqc . -c ${multiqc_config}
        """
}
