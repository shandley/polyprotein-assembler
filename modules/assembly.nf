process ASSEMBLY {
    tag "${sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_contigs.fa"), emit: contigs
    path "${sample_id}_assembly_log.txt"
    
    script:
    if (params.assembly_tool == 'megahit')
        """
        megahit \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            -o ${sample_id}_megahit \
            -t ${params.threads} \
            --min-contig-len ${params.min_contig_length} \
            --out-prefix ${sample_id}
            
        mv ${sample_id}_megahit/${sample_id}.contigs.fa ${sample_id}_contigs.fa
        cp ${sample_id}_megahit/log ${sample_id}_assembly_log.txt
        """
    else if (params.assembly_tool == 'spades')
        """
        spades.py \
            -1 ${reads[0]} \
            -2 ${reads[1]} \
            -o ${sample_id}_spades \
            -t ${params.threads} \
            --only-assembler
            
        # Filter contigs by length
        awk -v min_len=${params.min_contig_length} 'BEGIN {RS=">" ; ORS=""} length(\$2) >= min_len {print ">"\\$0}' \
            ${sample_id}_spades/contigs.fasta > ${sample_id}_contigs.fa
            
        cp ${sample_id}_spades/spades.log ${sample_id}_assembly_log.txt
        """
    else
        error "Invalid assembly tool: ${params.assembly_tool}. Choose either 'megahit' or 'spades'."
}