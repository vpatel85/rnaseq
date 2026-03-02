#!/usr/bin/env nextflow

/*
========================================================================================
    Basic RNA-seq Pipeline
========================================================================================
    A simple RNA-seq workflow demonstrating core analysis steps
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETERS
========================================================================================
*/

params.input          = null
params.genome_fasta   = null
params.genome_gtf     = null
params.outdir         = './results'
params.star_index     = null
params.skip_trimming  = false
params.help           = false

// Show help message
if (params.help) {
    log.info """
    Basic RNA-seq Pipeline
    ======================
    
    Usage:
        nextflow run main.nf --input samplesheet.csv --genome_fasta genome.fa --genome_gtf genes.gtf
    
    Required Arguments:
        --input          Path to samplesheet (CSV format: sample,fastq_1,fastq_2)
        --genome_fasta   Path to genome FASTA file
        --genome_gtf     Path to genome GTF annotation file
    
    Optional Arguments:
        --outdir         Output directory (default: ./results)
        --star_index     Path to prebuilt STAR index (if available)
        --skip_trimming  Skip adapter trimming step
    """.stripIndent()
    exit 0
}

// Validate required parameters
if (!params.input) {
    error "Please provide --input samplesheet"
}
if (!params.genome_fasta) {
    error "Please provide --genome_fasta"
}
if (!params.genome_gtf) {
    error "Please provide --genome_gtf"
}

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process FASTQC {
    tag "$meta.id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    
    script:
    def prefix = meta.id
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

process TRIMGALORE {
    tag "$meta.id"
    publishDir "${params.outdir}/trimgalore", mode: 'copy'
    
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_val_{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("*_trimming_report.txt"), emit: log
    
    script:
    def prefix = meta.id
    """
    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --gzip \\
        ${reads[0]} ${reads[1]}
    """
}

process STAR_INDEX {
    publishDir "${params.outdir}/star_index", mode: 'copy'
    
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_0'
    
    input:
    path fasta
    path gtf
    
    output:
    path "star_index", emit: index
    
    script:
    """
    mkdir star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100
    """
}

process STAR_ALIGN {
    tag "$meta.id"
    publishDir "${params.outdir}/star", mode: 'copy'
    
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_0'
    
    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    
    output:
    tuple val(meta), path("${meta.id}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}_Log.final.out"), emit: log
    tuple val(meta), path("${meta.id}_ReadsPerGene.out.tab"), emit: counts
    
    script:
    def prefix = meta.id
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${prefix}_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM NM \\
        --quantMode GeneCounts
    """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    publishDir "${params.outdir}/star", mode: 'copy'
    
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam_bai
    
    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

process FEATURECOUNTS {
    tag "$meta.id"
    publishDir "${params.outdir}/featurecounts", mode: 'copy'
    
    container 'quay.io/biocontainers/subread:2.0.6--he4a0461_2'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    
    output:
    tuple val(meta), path("${meta.id}_gene_counts.txt"), emit: counts
    tuple val(meta), path("${meta.id}_gene_counts.txt.summary"), emit: summary
    
    script:
    def prefix = meta.id
    """
    featureCounts \\
        -p \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o ${prefix}_gene_counts.txt \\
        ${bam}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    container 'quay.io/biocontainers/multiqc:1.25.2--pyhdfd78af_0'
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"
    
    script:
    """
    multiqc .
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    
    // Parse input samplesheet
    ch_input = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            [meta, reads]
        }
    
    // Run FastQC on raw reads
    FASTQC(ch_input)
    
    // Adapter trimming (optional)
    if (!params.skip_trimming) {
        TRIMGALORE(ch_input)
        ch_trimmed = TRIMGALORE.out.reads
    } else {
        ch_trimmed = ch_input
    }
    
    // Build or use existing STAR index
    if (params.star_index) {
        ch_star_index = channel.fromPath(params.star_index)
    } else {
        STAR_INDEX(
            channel.fromPath(params.genome_fasta),
            channel.fromPath(params.genome_gtf)
        )
        ch_star_index = STAR_INDEX.out.index
    }
    
    // Align reads with STAR
    STAR_ALIGN(
        ch_trimmed,
        ch_star_index,
        channel.fromPath(params.genome_gtf)
    )
    
    // Index BAM files
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
    
    // Quantify gene expression with featureCounts
    FEATURECOUNTS(
        SAMTOOLS_INDEX.out.bam_bai,
        channel.fromPath(params.genome_gtf)
    )
    
    // Collect all QC outputs for MultiQC
    ch_multiqc = channel.empty()
        .mix(FASTQC.out.zip.map { meta, files -> files })
        .mix(STAR_ALIGN.out.log.map { meta, files -> files })
        .mix(FEATURECOUNTS.out.summary.map { meta, files -> files })
    
    if (!params.skip_trimming) {
        ch_multiqc = ch_multiqc.mix(TRIMGALORE.out.log.map { meta, files -> files })
    }
    
    MULTIQC(ch_multiqc.collect())
}

workflow.onComplete {
    log.info """
    Pipeline completed!
    Status:    ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:  ${workflow.duration}
    Output:    ${params.outdir}
    """
}
