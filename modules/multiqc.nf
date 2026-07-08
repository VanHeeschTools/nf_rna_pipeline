process MULTIQC {
    label "multiqc"

    input:
        path multiqc_files, stageAs: "?/*"
        path multiqc_config

    output:
        path "*multiqc_report.html", emit: multiqc_report

    script:
        """
        multiqc . -c ${multiqc_config}
        """
}
