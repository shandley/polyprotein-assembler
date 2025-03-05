process FASTQC {
    tag "${sample_id}"
    label 'process_low'
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "fastqc_${sample_id}_logs"
    
    script:
    """
    mkdir fastqc_${sample_id}_logs
    
    fastqc -o fastqc_${sample_id}_logs \
        -t ${params.threads} \
        -q ${reads}
    """
}