# Nextflow RNA-seq Pipeline

This pipeline is designed for the analysis of RNA-seq data, covering various modules such as alignment, assembly, QC, fusion detection, and expression analysis.

## Index
- [Overview](#overview)
- [Usage](#usage)
  - [Default behaviour and run modes](#default-behaviour-and-run-modes)
- [Inputs](#inputs)
  - [Samplesheet](#samplesheet)
  - [Parameter specification](#parameter-specification)
  - [Optional inputs for modular execution](#optional-inputs-for-modular-execution)
- [Outputs](#outputs)
- [Support and Contributions](#support-and-contributions)

## Overview

1. **Quality Control and Adapter Trimming**: Utilizes FastP for assessing the quality of raw sequencing data and trims adapters and filters low-quality reads. Additionaly confirms strandedness for downstream transcriptome assembly.
2. **Alignment**: Aligns processed reads to a reference genome using STAR.
3. **Transcriptome Assembly**: Assembles transcripts of individual samples with stringtie, merges all samples and checks against a reference with GFFcompare, and finally filters the merged transcriptome.
4. **Fusions**: Detects gene fusions using Arriba, optionally takes vcf files from WGS to refine confidence thresholds.
5. **Expression Analysis**: Quantifies gene expression using Salmon.

## Usage

1. **Configure Parameters**:

    Customize the parameters using the `params.config` file to suit your analysis.

    - **Inputs**: Parameters are specified in `params.config`. See an example [here](https://github.com/mars13/nf_rna_pipeline/blob/dev/test/documentation/test_params.config). Alternatively, they can be provided as command-line arguments if double-dashed (`--samplesheet`). A detailed description of the required inputs can be found in the [Inputs](#inputs) section.

    - **Modules**: Adjust settings in `params.config` to toggle pipeline modules on and off (`qc`, `align`, `assembly`, `merge`, `fusion`, `expression`). All steps are set to `true` by default.
    
    ```
        params {
            //Toggle pipeline steps
            qc = false
            align = false
            assembly = false
            merge = true
        }
    ```

    - **Reference and other files**: Define paths to reference files, optional inputs, and other settings in `params.config`.

    - **Containers**: Modify process containers in `params.config` for tools like Trimmomatic, STAR, StringTie 2/3, etc., as per your environment. Note: This is expected to be provided with the pipeline in the future.

2. **Pull repository**:
    To keep up with the latest version pull the repository if it has already been downloaded.

    ```bash
    nextflow pull mars13/nf_rna_pipeline
    ```

3. **Run the Pipeline**:

    ```bash
    nextflow run mars13/nf_rna_pipeline -c params.config -profile [local/slurm]
    ```

    Additional nextflow run options can be provided. See [nextflow docs](https://www.nextflow.io/docs/latest/cli.html#run) for more information.

### Default behaviour and run modes

By default all module toggles are all set to true. 

The pipeline supports mode presets (`params.run_mode` or `--run_mode` from CLI) that change module defaults when the user did not explicitly set the individual toggles.

- qc
  - Intended for running only the QC step.
  - Behaviour (applied only when the individual toggle was NOT set by the user):
    - qc: true
    - align: false
    - assembly: false
    - merge: false
    - fusions: false
    - expression: false
  - Use case: quickly check data quality and strandedness without running downstream analyses.

- qc_restart
  - Intended to skip QC and run downstream steps (useful after QC already completed).
  - Behaviour: 
    - qc: false
    - align, assembly, merge, fusions, expressions: as defined by the user or true by default

  - Use case: re-run alignment/assembly/expression/fusion modules while avoiding re-running QC.

- manual (default)
    - No module is automatically enabled or disabled; module execution is controlled entirely by individual toggle parameters (e.g., `--qc`, `--align`, `--assembly`, etc.).

How to use:
- Set mode in params.config or via CLI:

    ```
    nextflow run . --run_mode qc
    nextflow run . --run_mode qc_restart
    ```

- You can still override any module individually:

    ```
    nextflow run . --run_mode qc --fusions false
    ```

## Inputs

### Samplesheet
- The samplesheet is a critical input for the RNA pipeline. It is a structured CSV file that lists all samples and files to be processed. Each row represents a sample with detailed information required by the pipeline. Below is a breakdown of the expected columns and their constraints:

| **Column Name**  | **Description**                                                                     | **Required** | **Constraints**                                                                 |
|------------------|---------------------------------------------------------------------------------|--------------|---------------------------------------------------------------------------------|
| `subject_id`     | Unique identifier for the subject (e.g., patient).                              | Yes          | Must be a string without spaces.                                                |
| `sample_id`      | Unique identifier for the sample (e.g., sample barcode).                        | Yes          | Must be a string without spaces.                                                |
| `group_id`       | (Optional) Identifier for the sample group (e.g., cohort).                      | No           | Must be a string without spaces. If left empty, it must be omitted entirely.    |
| `sample_type`    | Type of sample: either `tumor` or `normal`.                                     | Yes          | Must be either `tumor` or `normal`.                                             |
| `sequence_type`  | Specifies the type of sequencing data: `rna` or `dna`.                          | Yes          | Must be either `rna` or `dna`.                                                  |
| `file_type`      | Format of the input files, e.g., `fastq`.   | Yes          | Must be one of the supported formats: `fastq`, `bam`, `cram`, `vcf`, `csv`, etc. *|
| `filename_1`     | Path to the first file (e.g., R1 FASTQ file for paired-end or single-end data).  | Yes          | Must be a valid file path with no spaces. File extension must match `file_type`. |
| `filename_2`     | (Optional) Path to the second file (e.g., R2 FASTQ file for paired-end data).    | No           | Must be a valid file path with no spaces. If not applicable, leave empty.       |
| `conversion_reference` | (Optional) Reference fasta used to create a CRAM file (`file_type` `cram`). Only used to convert the CRAM to fastq. | No | Must be an existing fasta with a `.fai` index. If empty, the reference is matched from `conversion_references`. |

(*) Note: Not all suported formats can be used within the pipeline. RNA input can be `fastq` or `cram` (see [CRAM input](#cram-input)); bams are still not implemented as input.

**Example samplesheet**:
```
subject_id,sample_id,group_id,sample_type,sequence_type,file_type,filename_1,filename_2
subject1,sampleA,cohort1,tumor,rna,fastq,/path/to/sampleA_R1.fastq.gz,/path/to/sampleA_R2.fastq.gz
subject2,sampleB,,normal,rna,bam,/path/to/sampleB.bam,
subject3,sampleC,cohort2,tumor,dna,vcf,/path/to/sampleC.vcf,
```

Notes:
- All file paths (filename_1 and filename_2) must not contain spaces and should have extensions that match the declared file_type.
- Subject and sample ids must contain charaters (only numeric values will not be read properly).
- For single-end data or non-FASTQ files (e.g., BAM, CRAM, VCF), filename_2 can be omitted or left empty. For CRAM rows filename_2 is ignored: no `.crai` index is needed, because the whole file is converted.

### CRAM input

RNA-seq CRAM files (e.g. from a biobank) can be given directly in the samplesheet with `file_type` set to `cram`. Before QC, each CRAM is converted to paired fastq files with `samtools collate | samtools fastq`, after which the pipeline runs as it does for fastq input. CRAM and fastq rows can be mixed in one samplesheet.

Decoding a CRAM needs the reference fasta it was created with (the conversion reference). This reference is **only** used to convert the CRAM to fastq. All analysis steps after the conversion (alignment, assembly, expression, fusions) use the pipeline references (`reference_gtf`, `star_index_basedir`, `masked_fasta`, ...), so CRAM samples from other groups are analysed with the same references as fastq samples. The conversion reference is selected as follows:

1. If the `conversion_reference` column is filled for the sample, that fasta is used.
2. Otherwise the fastas listed in `conversion_references` are tried. The fasta whose sequence names and M5 checksums match all `@SQ` lines of the CRAM header is used, so it doesn't matter where the CRAM was originally created or whether that path is accessible.
3. If no fasta matches, the pipeline stops and reports the reference path stored in the CRAM header.

A reference given in the `conversion_reference` column is checked against the CRAM header in the same way. Every fasta needs a `.fai` index (`samtools faidx`). The selected reference of each sample is printed in the log (`cram input : <sample> uses conversion reference <fasta>`).

Example `conversion_references` file:
```
/path/to/references/GRCh38_v31_CTAT/ref_genome.fa
/path/to/references/GRCh38_v29_CTAT/ref_genome.fa
```

Notes:
- Only paired-end CRAM files are supported. CRAM rows count as paired-end in the paired/single-end check, so they can't be mixed with single-end fastq rows.
- Reads without a mate (singletons) are not used.
- Converted fastq files are only kept in the work directory unless `save_converted_fastq = true`, in which case they are published to `{outdir}/cram_to_fastq/`.
- Conversion runs in the `cramHeaderSq`, `referenceSq` and `cramToFastq` processes, which use the `samtools` label for the container.
- A stub test of the CRAM input handling (no containers needed) can be run from the repository root: `nextflow run . -stub -c test/cram/stub_test.config`.

### Parameter specification
| **Parameter**               | **Description**                                                                                     | **Type**   | **Required** | **Used by Module**                | **Default**                          |
|-----------------------------|-----------------------------------------------------------------------------------------------------|------------|--------------|------------------------------------|--------------------------------------|
| **Path specifications**      |                                                                                                     |            |              |                                    |                                      |
| `input`                     | Path to the samplesheet (CSV) file containing sample information.                                    | string     | Yes          |                                    |                                      |
| `project_folder`            | Path to the folder where outputs and logs will be stored.                                            | string     | No           | All modules                        | `.`                                  |
| `outdir`                    | Name of the output folder.                                                                          | string     | No           | All modules                        | `${project_folder}/output`           |
| **Reference files**          |                                                                                                     |            |              |                                    |                                      |
| `reference_gtf`             | Path to the reference annotation GTF file.                                                          | string     | Yes          |                                    |                                      |
| `reference_transcriptome`    | Path to the reference transcriptome assembly.                                                       | string     | No           | `align`, `expression`, `assembly`  |                                      |
| `star_index_basedir`        | Path to the STAR index folder.                                                                       | string     | No           | `align`                            |                                      |
| `refseq_gtf`                | Path to the RefSeq GTF file.                                                                        | string     | No           | `align`, `expression`              |                                      |
| `masked_fasta`              | Path to the masked FASTA file.                                                                      | string     | No           | `align`                            |                                      |
| `species`                   | Species name.                                                                                       | string     | No           | All modules                        | `Homo_sapiens`                       |
| `genome_version`            | Genome version.                                                                                     | string     | No           | All modules                        | `GRCh38`                             |
| `annot_version`             | Annotation version.                                                                                 | string     | No           | `build_annotation`                 | `102`                                |
| `twobit`                    | Path to the 2bit file.                                                                              | string     | No           | `assembly`                         |                                      |
| `kallisto_index`            | Path to the Kallisto index file.                                                                    | string     | No           | `qc`                               |                                      |
| `arriba_reference`          | Path to the pre-built Arriba reference folder.                                                      | string     | No           | `fusions`                          |                                      |
| **Modules**                  |                                                                                                     |            |              |                                    |                                      |
| `run_mode`                  | Pipeline mode preset that affects module defaults ("manual", "qc", "qc_restart").              | string     | No           | All modules                        | `manual`                               |
| `qc`                        | Toggle for the quality control module.                                                              | boolean    | No           | `qc`                               | `true`                               |
| `align`                     | Toggle for the alignment module.                                                                    | boolean    | No           | `align`                            | `true`                               |
| `assembly`                  | Toggle for the assembly module.                                                                     | boolean    | No           | `assembly`                         | `true`                               |
| `expression`                | Toggle for the expression quantification module.                                                    | boolean    | No           | `expression`                       | `true`                               |
| `fusions`                   | Toggle for the fusion detection module.                                                             | boolean    | No           | `fusions`                          | `true`                               |
| `merge`                     | Toggle for merging individual GTF files into a single file.                                         | boolean    | No           | `merge`                             | `true`                               |
| `build_annotation`          | Toggle for building a custom annotation from merged GTFs.                                           | boolean    | No           | `build_annotation`                 | `true`                               |
| **Extended pipeline options**|                                                                                                     |            |              |                                    |                                      |
| `custom_annotation`         | Formats for custom annotation outputs (e.g., "orfquant").                                           | string     | No           | `build_annotation`                 | `"orfquant"`                         |
| `expression_mode`           | Expression quantification mode for tools like Salmon quant (sq) and FeatureCounts (fc). Can be combined (e.g., "sqfc").              | string     | No           | `expression`                       | `"sqfc"`                             |
| `strandedness_check`        | Number of reads to use for strandedness check.                           | integer     | No           | `assembly`                         | 1000000      |
| `chr_exclusion_list`        | Path to the list of chromosomes to exclude during transcriptome assembly.                           | string     | No           | `assembly`                         | `assets/chr_exclusion_list.txt`      |
| `output_basename`           | Basename for merged GTF, custom annotation, and expression files.                                   | string     | No           | `merge`, `build_annotation`        | `transcriptome_full`                 |
| **Optional path definitions**|                                                                                                     |            |              |                                    |                                      |
| `resource_folder`           | Path to the annotation resource folder.                                                             | string     | No           | All modules                        |                                      |
| `container_folder`          | Path to the container location.                                                                     | string     | No           | All modules                        |                                      |
| `sample_gtf_list`           | Path to a file containing a list of individual sample GTFs for merging.                             | string     | No           | `merge`                            |                                      |
| `strand_info`               | Strandedness information, if available.                                                             | string     | No           | `assembly`                         |                                      |
| `conversion_references`     | Text file listing candidate reference fastas used to convert CRAM input to fastq, one path per line (see [CRAM input](#cram-input)). Not used for the analysis itself. | string | No | CRAM input |  |
| `save_converted_fastq`      | Publish fastq files converted from CRAM input to `{outdir}/cram_to_fastq/`.                        | boolean    | No           | CRAM input                         | `false`                              |
| **Maximum specifications**   |                                                                                                     |            |              |                                    |                                      |
| `max_memory`                | Maximum memory allocation per process.                                                              | string     | No           | All modules                        |                                      |
| `max_cpus`                  | Maximum CPU allocation per process.                                                                 | integer    | No           | All modules                        |                                      |
| `max_time`                  | Maximum time allowed per process.                                                                   | string     | No           | All modules                        |                                      |                      |

### Optional inputs for modular execution

- `strand_info`: When not running the QC step, the path to the file containing strand information for each sample can be provided. The file is expected to be plain text file (tab separated, no headers) with two columns: sample_id and strand information (`rf` for first strand,`fr` for second strand or `unstranded`). If `strand_info` is not provided the software will check the default location: `"${outdir}/check_strandedness/strandedness_all.txt"`. Example: 
```
sample_01 rf
sample_02 unstranded
```
- `bam_files`: When not running the ALIGN step before ASSEMBLY, the path to bam files can be provided as a regular expression. If `bam_files` is not provided the software will check the default location: `"${outdir}/star/**/*.Aligned.sortedByCoord.out.bam"`.
- `sample_gtf_list`: For merging GTF files from individual samples without prior assembly steps, you are require to provide a file containing the paths to individual GTFs. Example:

```
./analysis/stringtie/sample_01/sample_01.gtf
./analysis/stringtie/sample_02/sample_02.gtf
```


## Outputs

The pipeline generates output files including quality reports, trimmed reads, alignment results and assembled transcripts. The output directory is specified by the parameter `outdir` and has the following structure:

```
{outdir}
├── check_strandedness/
├── customannotation/
├── fastp/
├── gffcompare/
├── salmon/
├── star/
└── stringtie/
```

Additionally it produces Nextflow execution reports in `{project_folder}/log`

## Support and Contributions

- For issues or questions, [create an issue](https://github.com/mars13/nf_rna_pipeline/issues).
- Contributions to enhance this pipeline are welcomed through pull requests.
