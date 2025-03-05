process TRIMMING {
    tag "${sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("trimmed/${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: trimmed_reads
    path "trimmed/${sample_id}_trim.log"
    
    script:
    """
    mkdir -p trimmed
    
    trimmomatic PE \
        -threads ${params.threads} \
        -phred33 \
        ${reads[0]} ${reads[1]} \
        trimmed/${sample_id}_R1_trimmed.fastq.gz trimmed/${sample_id}_R1_unpaired.fastq.gz \
        trimmed/${sample_id}_R2_trimmed.fastq.gz trimmed/${sample_id}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${params.adapter_file}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        2> trimmed/${sample_id}_trim.log
    """
}