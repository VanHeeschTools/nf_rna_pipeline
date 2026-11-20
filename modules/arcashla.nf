// HLA typing with arcasHLA (https://github.com/RabadanLab/arcasHLA)
// arcasHLA names its outputs after the part of the file name before the first dot,
// so the inputs are renamed to the sample id (sample ids with a dot are rejected in checkInputFiles)

// Define process to check that the IMGT/HLA database is available to arcasHLA, before any sample is processed.
// arcasHLA reads the database from the dat/ folder of its installation and tries to download it there when
// missing, which fails in a read-only container. A prepared dat/ folder can be bind-mounted instead.
process arcasReference {
    label "arcashla"

    output:
        path "arcashla_versions.txt", emit: versions // arcasHLA installation and IMGT/HLA database version

    script:
        """
        # Locate the arcasHLA installation (the folder with scripts/ and dat/).
        # Bioconda installs a wrapper script in bin/ that calls <share>/arcas-hla-<version>/arcasHLA.
        arcas_bin=\$(command -v arcasHLA)
        arcas_main=\$(grep -m1 -oE '/[^ ]+/arcasHLA( |\$)' "\${arcas_bin}" | tr -d ' ' || true)
        arcas_dir=\$(dirname "\$(readlink -f "\${arcas_main:-\${arcas_bin}}")")
        if [ ! -f "\${arcas_dir}/scripts/reference.py" ]; then
            echo "ERROR: arcasHLA installation not found (checked \${arcas_dir})" >&2
            exit 1
        fi

        if [ ! -s "\${arcas_dir}/dat/IMGTHLA/hla.dat" ] || [ ! -d "\${arcas_dir}/dat/ref" ]; then
            echo "ERROR: no IMGT/HLA database found in \${arcas_dir}/dat of the arcasHLA container." >&2
            echo "Use an arcasHLA image with a built-in database (containers/arcashla/arcashla.def)," >&2
            echo "or bind-mount a prepared dat folder to \${arcas_dir}/dat, e.g. with" >&2
            echo "  containerOptions = '-B /path/to/arcashla_dat:\${arcas_dir}/dat' (see README, HLA typing)" >&2
            exit 1
        fi

        # The bioconda package contains Git LFS placeholders instead of the fasta files used by arcasHLA customize
        if grep -l "git-lfs.github.com" "\${arcas_dir}"/dat/ref/*.fasta >&2 2>/dev/null; then
            echo "ERROR: the fasta files above in \${arcas_dir}/dat/ref are Git LFS placeholders, not the real files." >&2
            echo "Use an image built from containers/arcashla/arcashla.def, which downloads them (see README, HLA typing)." >&2
            exit 1
        fi

        # IMGT/HLA release from release_version.txt ("# version: IPD-IMGT/HLA 3.xx.0"), else the git commit
        imgt_version=\$(awk -F': *' 'tolower(\$1) ~ /^# *version\$/ {print \$2; exit}' "\${arcas_dir}/dat/IMGTHLA/release_version.txt" 2>/dev/null || true)
        if [ -z "\${imgt_version}" ]; then
            imgt_version=\$(git -C "\${arcas_dir}/dat/IMGTHLA" rev-parse HEAD 2>/dev/null || echo "unknown")
        fi
        printf 'arcasHLA installation\\t%s\\nIMGT/HLA database\\t%s\\n' "\${arcas_dir}" "\${imgt_version}" > arcashla_versions.txt
        """

    stub:
        """
        printf 'arcasHLA installation\\tstub\\nIMGT/HLA database\\tstub\\n' > arcashla_versions.txt
        """
}

// Define process to extract chromosome 6 and unmapped reads from a sorted bam
process arcasExtract {
    label "arcashla"
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(bam, stageAs: "input/*"), path(bai, stageAs: "input/*") // Tuple of sample id, sorted bam and bam index (empty list if not available)
        val paired_end                                                                    // Bool, true if data is paired end

    output:
        tuple val(sample_id), path("${sample_id}.extracted*.fq.gz"), emit: fastq // Tuple of sample id and extracted HLA reads

    script:
        def single = paired_end ? "" : "--single"
        def link_bai = bai ? "ln -s \"${bai}\" \"${sample_id}.bam.bai\"" : "# No bam index available, arcasHLA indexes the bam"
        """
        # Name the bam after the sample id, arcasHLA names its outputs after the bam
        ln -s "${bam}" "${sample_id}.bam"
        ${link_bai}
        mkdir -p tmp

        arcasHLA extract \
            "${sample_id}.bam" \
            --outdir . \
            --temp tmp/ \
            --threads ${task.cpus} \
            --unmapped \
            ${single} \
            --verbose

        rm -rf tmp
        """

    stub:
        def read_files = paired_end ? "1 2" : ""
        """
        for read in ${read_files ?: '""'}; do
            suffix=\${read:+.\${read}}
            printf '@read\\nACGT\\n+\\nIIII\\n' | gzip > "${sample_id}.extracted\${suffix}.fq.gz"
        done
        """
}

// Define process to predict the HLA genotype from the extracted reads
process arcasGenotype {
    label "arcashla"
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(fastq) // Tuple of sample id and extracted HLA reads
        val paired_end                    // Bool, true if data is paired end

    output:
        tuple val(sample_id), path("${sample_id}.genotype.json"), emit: genotype // Tuple of sample id and predicted HLA genotype
        path "${sample_id}.genes.json", emit: genes                               // Read counts per HLA gene, used for the cohort table
        path "${sample_id}.{genotype.log,alignment.p}"                            // Log and alignment for publishing

    script:
        def single = paired_end ? "" : "--single"
        def reads = [fastq].flatten().collect { read_file -> read_file.name }.sort().join(' ')
        """
        mkdir -p tmp

        arcasHLA genotype \
            ${reads} \
            --outdir . \
            --temp tmp/ \
            --threads ${task.cpus} \
            ${single} \
            --verbose

        rm -rf tmp
        """

    stub:
        """
        echo '{"A": ["A*01:01:01", "A*02:01:01"], "B": ["B*07:02:01", "B*08:01:01"]}' > "${sample_id}.genotype.json"
        echo '{"A": 100, "B": 80}' > "${sample_id}.genes.json"
        touch "${sample_id}.genotype.log" "${sample_id}.alignment.p"
        """
}

// Define process to quantify HLA allele expression using a sample specific HLA reference
process arcasQuant {
    label "arcashla"
    tag "${sample_id}"

    input:
        tuple val(sample_id), path(fastq), path(genotype) // Tuple of sample id, extracted HLA reads and predicted HLA genotype
        val paired_end                                    // Bool, true if data is paired end

    output:
        path "quant/${sample_id}.quant*", emit: quant                                                    // HLA allele and gene expression tables
        tuple val(sample_id), path("quant/${sample_id}.quant_skipped_genes.tsv"), emit: skipped // Genes that could not be quantified

    script:
        def single = paired_end ? "" : "--single"
        def reads = [fastq].flatten().collect { read_file -> read_file.name }.sort().join(' ')
        """
        mkdir -p tmp ref quant

        # Locate the arcasHLA installation (see arcasReference)
        arcas_bin=\$(command -v arcasHLA)
        arcas_main=\$(grep -m1 -oE '/[^ ]+/arcasHLA( |\$)' "\${arcas_bin}" | tr -d ' ' || true)
        arcas_dir=\$(dirname "\$(readlink -f "\${arcas_main:-\${arcas_bin}}")")

        # Only genes whose typed alleles all have a cDNA sequence in arcasHLA's fixed dat/ref/cDNA.json can be
        # quantified, other genes would make quant fail or be quantified incompletely (see arcashla_quant_genes.py)
        genes=\$(arcashla_quant_genes.py "${genotype}" "\${arcas_dir}" "quant/${sample_id}.quant_skipped_genes.tsv")
        if [ -z "\${genes}" ]; then
            echo "No HLA gene of ${sample_id} can be quantified, see quant/${sample_id}.quant_skipped_genes.tsv" >&2
            rm -rf tmp ref
            exit 0
        fi

        # Build HLA reference with the alleles of this sample
        arcasHLA customize \
            --genotype ${genotype} \
            --subject ${sample_id} \
            --genes "\${genes}" \
            --outdir ref \
            --temp tmp/ \
            --threads ${task.cpus} \
            --verbose

        # Quantify expression of the sample's HLA alleles
        arcasHLA quant \
            --ref ref/${sample_id} \
            --sample ${sample_id} \
            --outdir quant \
            --temp tmp/ \
            --threads ${task.cpus} \
            ${single} \
            --verbose \
            ${reads}

        rm -rf tmp
        """

    stub:
        """
        mkdir -p quant
        touch quant/${sample_id}.quant.tsv quant/${sample_id}.quant.alleles.tsv quant/${sample_id}.quant.genes.tsv
        printf 'gene\ttyped_alleles\tmissing_in_arcashla_cdna\treason\n' > quant/${sample_id}.quant_skipped_genes.tsv
        """
}

// Define process to merge the genotypes of all samples into cohort tables
process arcasMerge {
    label "arcashla"

    input:
        path json_files, stageAs: "json/*" // genotype.json and genes.json files of all samples
        val run_name                       // Prefix of the output tables

    output:
        path "${run_name}.genotypes.tsv", emit: genotypes // HLA genotype per sample
        path "${run_name}.genes.tsv", emit: genes         // Read counts per HLA gene per sample

    script:
        """
        arcasHLA merge \
            --indir json \
            --outdir . \
            --run ${run_name} \
            --verbose
        """

    stub:
        """
        touch ${run_name}.genotypes.tsv ${run_name}.genes.tsv
        """
}
