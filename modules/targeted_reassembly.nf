process TARGETED_REASSEMBLY {
    tag "${sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/targeted_reassembly", mode: 'copy'
    
    input:
    tuple val(sample_id), path(partial_polyproteins)
    tuple val(sample_id_reads), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_reassembled_contigs.fa"), emit: reassembled_contigs
    path "${sample_id}_reassembly_log.txt"
    
    script:
    """
    # Extract contig IDs containing partial polyproteins
    python3 << EOF
import json
import os
from Bio import SeqIO

# Load the partial polyprotein hits
with open("${partial_polyproteins}") as f:
    partial_hits = json.load(f)

# Write contig IDs to a file
with open("partial_contigs.txt", "w") as out:
    for contig_id in partial_hits.keys():
        out.write(f"{contig_id}\\n")
EOF

    # Create a script to extract reads mapping to the partial contigs
    cat << 'EOT' > extract_reads.sh
#!/bin/bash
set -e

# Map reads to original assembly
bwa index ${params.assembly_dir}/${sample_id}_contigs.fa
bwa mem -t ${params.threads} ${params.assembly_dir}/${sample_id}_contigs.fa ${reads[0]} ${reads[1]} | \\
    samtools view -b -F 4 - | \\
    samtools sort -o mapped.bam -

# Extract read IDs mapping to partial contigs
samtools index mapped.bam
while read contig; do
    samtools view mapped.bam "\$contig" | cut -f1 >> read_ids.txt
done < partial_contigs.txt

# Get unique read IDs
sort read_ids.txt | uniq > uniq_read_ids.txt

# Extract read sequences
seqtk subseq ${reads[0]} uniq_read_ids.txt > partial_R1.fastq
seqtk subseq ${reads[1]} uniq_read_ids.txt > partial_R2.fastq
EOT

    # Execute the extraction script
    chmod +x extract_reads.sh
    ./extract_reads.sh
    
    # Run assembly with smaller k-mers for improved sensitivity
    megahit \
        -1 partial_R1.fastq \
        -2 partial_R2.fastq \
        -o ${sample_id}_targeted \
        -t ${params.threads} \
        --k-min 21 \
        --k-max 101 \
        --k-step 10 \
        --min-contig-len ${params.min_contig_length} \
        --out-prefix ${sample_id}_targeted
        
    mv ${sample_id}_targeted/${sample_id}_targeted.contigs.fa ${sample_id}_reassembled_contigs.fa
    cp ${sample_id}_targeted/log ${sample_id}_reassembly_log.txt
    """
}