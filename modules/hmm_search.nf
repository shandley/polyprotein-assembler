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
    # Check if pyhmmer is available (better for Apple Silicon)
    if python3 -c "import pyhmmer" 2>/dev/null; then
        echo "Using pyhmmer for HMM search (optimized for ARM64/Apple Silicon)"
        # Run HMM search using pyhmmer
        python3 << EOF
import json
import os
import pyhmmer
from pyhmmer.plan7 import HMMFile
from Bio import SeqIO

# Load the HMM database
hmms = []
with pyhmmer.plan7.HMMFile("${hmm_db}") as hmm_file:
    for hmm in hmm_file:
        hmms.append(hmm)

# Load the protein sequences
sequences = []
for record in SeqIO.parse("${translated_contigs}", "fasta"):
    seq = pyhmmer.easel.TextSequence(
        name=str(record.id).encode(),
        sequence=str(record.seq).encode()
    )
    sequences.append(seq)

# Run the search
hits = []
for hmm in hmms:
    # Create a pipeline
    pipeline = pyhmmer.plan7.Pipeline()
    
    # Set e-value thresholds
    pipeline.E = ${params.hmm_evalue}
    pipeline.domE = ${params.hmm_dom_evalue}
    
    # Search with the hmm
    results = pipeline.search_hmm(hmm, sequences)
    
    # Process the results
    for hit in results:
        for domain in hit.domains:
            hit_entry = {
                "target_name": hit.name.decode(),
                "accession": "",
                "query_name": hmm.name.decode(),
                "query_accession": hmm.accession.decode() if hmm.accession else "",
                "evalue": hit.evalue,
                "score": hit.score,
                "bias": hit.bias,
                "domain_evalue": domain.evalue,
                "domain_score": domain.score,
                "domain_bias": domain.bias,
                "hmm_from": domain.hmm_from,
                "hmm_to": domain.hmm_to,
                "ali_from": domain.env_from,
                "ali_to": domain.env_to,
                "env_from": domain.env_from,
                "env_to": domain.env_to,
                "coverage": (domain.hmm_to - domain.hmm_from + 1) / hmm.length,
                "contig_id": hit.name.decode().split("_frame")[0]
            }
            hits.append(hit_entry)

# Create tabular output for compatibility
with open("${sample_id}_hmmsearch.tbl", "w") as f:
    f.write("# hmmsearch :: search profile(s) against a sequence database\\n")
    f.write("# HMMER 3.3.2 (via pyhmmer)\\n")
    f.write("# Query: Query HMM\\n")
    f.write("# Target: Target sequence database\\n")
    f.write("# Option: -E ${params.hmm_evalue} --domE ${params.hmm_dom_evalue}\\n")
    f.write("#\\n")
    f.write("# Domain annotation for each hit (and alignments):\\n")
    f.write("#\\n")
    f.write("# Scores for sequence family classification (score includes all domains):\\n")
    f.write("# Model      Description                          Seq            Domain  score    E-value  score    E-value  hmm_from   hmm_to    ali_from   ali_to    env_from   env_to     acc\\n")
    f.write("# --------   -----------                          -------------- ------- -----    ------- -----    ------- --------   --------   --------   --------   --------   --------   ----\\n")
    
    for hit in hits:
        f.write(f"{hit['target_name']}\\t{hit['accession']}\\t{hit['query_name']}\\t{hit['query_accession']}\\t")
        f.write(f"0\\t0\\t{hit['evalue']:.2e}\\t{hit['score']:.1f}\\t{hit['bias']:.1f}\\t0\\t0\\t")
        f.write(f"{hit['domain_evalue']:.2e}\\t{hit['domain_score']:.1f}\\t{hit['domain_bias']:.1f}\\t")
        f.write(f"{hit['hmm_from']}\\t{hit['hmm_to']}\\t{hit['ali_from']}\\t{hit['ali_to']}\\t")
        f.write(f"{hit['env_from']}\\t{hit['env_to']}\\t0.00\\n")

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
    else
        echo "Using standard hmmsearch (pyhmmer not available)"
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
    fi
    """
}