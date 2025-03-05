process HMM_SEARCH {
    tag "${sample_id}"
    label 'process_high'
    
    publishDir "${params.outdir}/hmm_search", mode: 'copy'
    
    input:
    tuple val(sample_id), path(translated_contigs)
    path hmm_db
    
    output:
    tuple val(sample_id), path("${sample_id}_hmm_hits.json"), emit: hmm_hits
    path "${sample_id}_hmmsearch.tbl"
    
    script:
    """
    # Run hmmsearch against the polyprotein profiles
    hmmsearch --cpu ${params.threads} \
        --domtblout ${sample_id}_hmmsearch.tbl \
        -E ${params.hmm_evalue} \
        --domE ${params.hmm_dom_evalue} \
        ${hmm_db} \
        ${translated_contigs}
    
    # Convert the tabular output to JSON for easier processing
    python3 << EOF
import json
import os

hits = []
with open("${sample_id}_hmmsearch.tbl") as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        if len(fields) < 22:
            continue
            
        hit = {
            "target_name": fields[0],
            "accession": fields[1],
            "query_name": fields[3],
            "query_accession": fields[4],
            "evalue": float(fields[6]),
            "score": float(fields[7]),
            "bias": float(fields[8]),
            "domain_evalue": float(fields[12]),
            "domain_score": float(fields[13]),
            "domain_bias": float(fields[14]),
            "hmm_from": int(fields[15]),
            "hmm_to": int(fields[16]),
            "ali_from": int(fields[17]),
            "ali_to": int(fields[18]),
            "env_from": int(fields[19]),
            "env_to": int(fields[20]),
            "coverage": (int(fields[16]) - int(fields[15]) + 1) / float(fields[5]),
            "contig_id": fields[0].split("_frame")[0]
        }
        hits.append(hit)

# Group by contig
contig_hits = {}
for hit in hits:
    contig_id = hit["contig_id"]
    if contig_id not in contig_hits:
        contig_hits[contig_id] = []
    contig_hits[contig_id].append(hit)

with open("${sample_id}_hmm_hits.json", 'w') as f:
    json.dump(contig_hits, f, indent=2)
EOF
    """
}