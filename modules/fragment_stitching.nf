process FRAGMENT_STITCHING {
    tag "${sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/fragment_stitching", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reassembled_contigs)
    tuple val(sample_id_full), path(full_polyproteins)
    
    output:
    tuple val(sample_id), path("${sample_id}_improved_polyproteins.json"), emit: improved_polyproteins
    path "${sample_id}_stitching_report.html"
    
    script:
    """
    # Translate the reassembled contigs
    transeq -sequence ${reassembled_contigs} -outseq ${sample_id}_reassembled.faa -frame=6 -clean
    
    # Run hmmsearch on the reassembled contigs
    hmmsearch --cpu ${params.threads} \
        --domtblout ${sample_id}_reassembled_hmmsearch.tbl \
        -E ${params.hmm_evalue} \
        --domE ${params.hmm_dom_evalue} \
        ${params.hmm_db} \
        ${sample_id}_reassembled.faa
    
    # Stitch fragments together to create improved polyproteins
    python3 << EOF
import json
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from collections import defaultdict

# Parse the reassembled contigs hits
hits = []
with open("${sample_id}_reassembled_hmmsearch.tbl") as f:
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

# Load original full polyproteins
with open("${full_polyproteins}") as f:
    full_hits = json.load(f)

# Load translated sequences
reassembled_seqs = {record.id: record for record in SeqIO.parse("${sample_id}_reassembled.faa", "fasta")}

# Group hits by HMM profile
profile_hits = defaultdict(list)
for hit in hits:
    profile_hits[hit["query_name"]].append(hit)

# Create graph of overlapping fragments
improved_polyproteins = {}

for profile, profile_hit_list in profile_hits.items():
    if len(profile_hit_list) <= 1:
        continue
        
    # Sort hits by HMM position
    profile_hit_list.sort(key=lambda x: x["hmm_from"])
    
    # Create a graph of overlapping fragments
    G = nx.Graph()
    
    # Add nodes for each hit
    for i, hit in enumerate(profile_hit_list):
        G.add_node(i, hit=hit)
    
    # Add edges between overlapping fragments
    for i in range(len(profile_hit_list)):
        hit_i = profile_hit_list[i]
        for j in range(i+1, len(profile_hit_list)):
            hit_j = profile_hit_list[j]
            
            # Check for overlap in HMM coordinates
            if hit_i["hmm_to"] >= hit_j["hmm_from"] and hit_i["hmm_from"] <= hit_j["hmm_to"]:
                # Calculate overlap score
                overlap_len = min(hit_i["hmm_to"], hit_j["hmm_to"]) - max(hit_i["hmm_from"], hit_j["hmm_from"])
                overlap_score = overlap_len / min(hit_i["hmm_to"] - hit_i["hmm_from"], hit_j["hmm_to"] - hit_j["hmm_from"])
                
                if overlap_score >= ${params.min_overlap}:
                    G.add_edge(i, j, weight=overlap_score)
    
    # Find connected components (potential polyprotein fragments)
    connected_components = list(nx.connected_components(G))
    
    for i, component in enumerate(connected_components):
        if len(component) >= ${params.min_fragments}:
            # Sort the component nodes by HMM position
            component_sorted = sorted(component, key=lambda x: profile_hit_list[x]["hmm_from"])
            
            # Create stitched polyprotein
            component_key = f"{profile}_{i+1}"
            
            # Store the stitched polyprotein info
            component_hits = [profile_hit_list[idx] for idx in component_sorted]
            
            # Calculate completeness based on HMM coverage
            hmm_positions = []
            for hit in component_hits:
                hmm_positions.extend(range(hit["hmm_from"], hit["hmm_to"] + 1))
            
            unique_positions = len(set(hmm_positions))
            hmm_length = int(component_hits[0]["query_accession"].split('.')[1])
            completeness = unique_positions / hmm_length
            
            improved_polyproteins[component_key] = {
                "profile": profile,
                "hits": component_hits,
                "completeness": completeness,
                "fragment_count": len(component_hits),
                "total_score": sum(hit["score"] for hit in component_hits)
            }

# Generate a report
report_data = []
for key, polyprotein in improved_polyproteins.items():
    report_data.append({
        "Polyprotein ID": key,
        "Profile": polyprotein["profile"],
        "Completeness": f"{polyprotein['completeness']:.2f}",
        "Fragment Count": polyprotein["fragment_count"],
        "Total Score": f"{polyprotein['total_score']:.2f}"
    })

# Sort by completeness
report_df = pd.DataFrame(report_data)
if not report_df.empty:
    report_df = report_df.sort_values("Completeness", ascending=False)

# Generate HTML report
html = f"""
<html>
<head>
    <title>Polyprotein Fragment Stitching Report - {sample_id}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
        }}
        th, td {{
            text-align: left;
            padding: 8px;
            border: 1px solid #ddd;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        h1, h2 {{
            color: #333;
        }}
        .summary {{
            margin-bottom: 20px;
        }}
    </style>
</head>
<body>
    <h1>Polyprotein Fragment Stitching Report</h1>
    <div class="summary">
        <h2>Summary</h2>
        <p>Sample ID: {sample_id}</p>
        <p>Total improved polyproteins: {len(improved_polyproteins)}</p>
    </div>
    
    <h2>Improved Polyproteins</h2>
    {report_df.to_html(index=False) if not report_df.empty else "<p>No improved polyproteins found.</p>"}
</body>
</html>
"""

with open("${sample_id}_stitching_report.html", "w") as f:
    f.write(html)

# Save the improved polyproteins
with open("${sample_id}_improved_polyproteins.json", "w") as f:
    json.dump(improved_polyproteins, f, indent=2)
EOF
    """
}