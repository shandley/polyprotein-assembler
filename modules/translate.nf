process TRANSLATE_CONTIGS {
    tag "${sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/translated_contigs", mode: 'copy'
    
    input:
    tuple val(sample_id), path(contigs)
    
    output:
    tuple val(sample_id), path("${sample_id}_translated.faa"), emit: translated_contigs
    
    script:
    """
    # Translate contigs in all 6 reading frames using transeq from EMBOSS
    transeq -sequence ${contigs} -outseq ${sample_id}_6frame.faa -frame=6 -clean
    
    # Filter translated sequences to remove those with stop codons
    python3 << EOF
import os
from Bio import SeqIO

with open("${sample_id}_6frame.faa") as input_handle, open("${sample_id}_translated.faa", "w") as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        # Skip sequences with internal stop codons or too short
        if "*" not in record.seq[:-1] and len(record.seq) >= ${params.min_protein_length}:
            SeqIO.write(record, output_handle, "fasta")
EOF
    """
}