#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "*.fastq.gz"
params.outdir = "results"
params.genome = "/media/dimbo/10T/TAL_LAB/Genomes/hisat_index/mm10.fa"
params.genome_index = "/media/dimbo/10T/TAL_LAB/Genomes/hisat_index/mm10"
params.gsize = "mm"  // Genome size for MACS2
params.merge_groups = [:]  // Will be populated from config file

process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'
    conda "bioconda::fastqc=0.11.9"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -q -t ${task.cpus} ${reads}
    """
}

process TRIMMOMATIC {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmomatic", mode: 'copy'
    conda "bioconda::trimmomatic=0.39"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads

    script:
    """
    TrimmomaticSE -threads ${task.cpus} -phred33 \
        ${reads} \
        ${sample_id}_trimmed.fastq.gz \
        SLIDINGWINDOW:4:18 \
        LEADING:28 \
        TRAILING:28 
    """
}

process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc/trimmed", mode: 'copy'
    conda "bioconda::fastqc=0.11.9"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -q -t ${task.cpus} ${reads}
    """
}

process HISAT2_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/hisat", mode: 'copy'
    conda "bioconda::hisat2=2.2.1 bioconda::samtools=1.15"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    """
    hisat2 -p ${task.cpus} --no-spliced-alignment -x ${params.genome_index} -U ${reads} | \
    samtools view -@ ${task.cpus} -bS - | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    """
}

process SAMTOOLS_RMDUP {
    tag "$sample_id"
    publishDir "${params.outdir}/rmdup", mode: 'copy'
    conda "bioconda::samtools=1.15"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.rmdup.bam"), emit: rmdup_bam

    script:
    """
    samtools rmdup ${bam} ${sample_id}.rmdup.bam
    """
}

process SAMTOOLS_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered_bam", mode: 'copy'
    conda "bioconda::samtools=1.15"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), emit: filtered_bam

    script:
    """
    samtools view -@ ${task.cpus} -bh -q 30 -F 1024 ${bam} > ${sample_id}.filtered.bam
    """
}

process MERGE_BAM {
    tag "$group_id"
    publishDir "${params.outdir}/merged", mode: 'copy'
    conda "bioconda::samtools=1.15"

    input:
    tuple val(group_id), path(bams)

    output:
    tuple val(group_id), path("${group_id}.bam"), emit: merged_bam

    script:
    """
    samtools merge -@ ${task.cpus} ${group_id}.bam ${bams}
    """
}

process COUNT_READS {
    tag "$sample_id"
    conda "bioconda::samtools=1.15"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), stdout, emit: counts

    script:
    """
    samtools view -c ${bam} | tr -d '\n'
    """
}

process DOWNSAMPLE_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/downsampled", mode: 'copy'
    conda "bioconda::samtools=1.15"

    input:
    tuple val(sample_id), path(bam)
    val min_reads

    output:
    tuple val(sample_id), path("${sample_id}_downsampled.bam"), emit: downsampled_bam

    script:
    """
    total_reads=\$(samtools view -c ${bam})
    fraction=\$(echo "scale=4; ${min_reads} / \$total_reads" | bc)
    samtools view -@ ${task.cpus} -bs \$fraction ${bam} > ${sample_id}_downsampled.bam
    """
}

process SAMTOOLS_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/indexed", mode: 'copy'
    conda "bioconda::samtools=1.15"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

process BAM_TO_BIGWIG {
    tag "$sample_id"
    publishDir "${params.outdir}/bigwig", mode: 'copy'
    conda "bioconda::deeptools=3.5.1"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.bw"), emit: bigwig

    script:
    """
    bamCoverage -b ${bam} -o ${sample_id}.bw -p ${task.cpus}
    """
}

process MACS2_CALLPEAK {
    tag "$sample_id"
    publishDir "${params.outdir}/macs2", mode: 'copy'
    conda "bioconda::macs2=2.2.7.1"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_peaks_peaks.narrowPeak"), emit: peaks
    path "${sample_id}_peaks_*"

    script:
    """
    macs2 callpeak --treatment ${bam} \
                   --format BAM \
                   --gsize ${params.gsize} \
                   -n ${sample_id}_peaks \
                   -p 0.0000000000001

    # Ensure the output file exists with the expected name
    if [ -f "${sample_id}_peaks_peaks.narrowPeak" ]; then
        echo "Peak file found with expected name"
    elif [ -f "${sample_id}_peaks.narrowPeak" ]; then
        mv "${sample_id}_peaks.narrowPeak" "${sample_id}_peaks_peaks.narrowPeak"
        echo "Peak file renamed to match expected output"
    else
        echo "Peak file not found" >&2
        exit 1
    fi
    """

}

workflow {
    Channel
        .fromPath(params.reads)
        .map { file -> tuple(file.simpleName, file) }
        .set { read_ch }

    FASTQC_RAW(read_ch)
    trimmed_reads = TRIMMOMATIC(read_ch)
    FASTQC_TRIMMED(trimmed_reads)
    
    aligned_reads = HISAT2_ALIGN(trimmed_reads)
    rmdup_reads = SAMTOOLS_RMDUP(aligned_reads)
    filtered_reads = SAMTOOLS_FILTER(rmdup_reads)
    
    // Prepare channel for merging
    reads_for_merge = filtered_reads
        .map { sample_id, bam -> 
            def group = params.merge_groups.find { group, samples -> samples.contains(sample_id) }
            group ? tuple(group.key, bam) : tuple(sample_id, bam)
        }
        .groupTuple()

    // Merge BAM files
    merged_reads = MERGE_BAM(reads_for_merge)

    // Continue analysis only with merged files
    
    // Count reads in merged samples
    read_counts = COUNT_READS(merged_reads)
    
    // Find the minimum read count across merged samples
    min_reads = read_counts
        .map { it[1].toLong() }
        .min()
    
    // Downsample merged samples to the minimum read count
    downsampled_reads = DOWNSAMPLE_BAM(merged_reads, min_reads)

    // Index the downsampled BAM files
    indexed_reads = SAMTOOLS_INDEX(downsampled_reads)

    // Convert downsampled BAM to BigWig
    BAM_TO_BIGWIG(indexed_reads)

    // Call peaks using MACS2
    MACS2_CALLPEAK(downsampled_reads)
}
