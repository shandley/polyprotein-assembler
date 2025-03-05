process ANNOTATION {
    tag "${sample_id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/annotation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(scored_polyproteins)
    
    output:
    path "${sample_id}_annotated_polyproteins.fa"
    path "${sample_id}_annotations.gff"
    path "${sample_id}_annotation_report.html"
    
    script:
    """
    # Annotate the polyproteins with functional information
    python3 << EOF
import json
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from collections import defaultdict

# Load the scored polyproteins
with open("${scored_polyproteins}") as f:
    polyproteins = json.load(f)

# Create a FASTA file with the polyprotein sequences
polyprotein_records = []
annotations = []

for key, polyprotein in polyproteins.items():
    # Extract hit information
    profile = polyprotein["profile"]
    completeness = polyprotein["completeness"]
    completeness_class = polyprotein["completeness_class"]
    
    # Load the original sequence for each hit
    hits = polyprotein["hits"]
    
    for i, hit in enumerate(hits):
        target_name = hit["target_name"]
        contig_id = hit["contig_id"]
        
        # Create a description with annotation info
        description = f"profile={profile} completeness={completeness:.2f} class={completeness_class}"
        
        # Create a new record ID
        record_id = f"{key}_fragment_{i+1}"
        
        # Note: In a real implementation, we would extract the actual sequence
        # from the original contigs. For this example, we'll create a placeholder.
        seq = Seq("A" * 100)  # Placeholder sequence
        record = SeqRecord(seq, id=record_id, description=description)
        polyprotein_records.append(record)
        
        # Create annotation entries
        annotations.append({
            "seqid": record_id,
            "source": "PolyproteinAssembler",
            "type": "CDS",
            "start": 1,
            "end": 100,  # Placeholder
            "score": hit["score"],
            "strand": "+",
            "phase": 0,
            "attributes": f"ID={record_id};Profile={profile};Completeness={completeness:.2f};Class={completeness_class}"
        })

# Write the polyprotein sequences to a FASTA file
with open("${sample_id}_annotated_polyproteins.fa", "w") as f:
    SeqIO.write(polyprotein_records, f, "fasta")

# Write the annotations to a GFF file
with open("${sample_id}_annotations.gff", "w") as f:
    f.write("##gff-version 3\\n")
    for annotation in annotations:
        f.write(f"{annotation['seqid']}\\t{annotation['source']}\\t{annotation['type']}\\t")
        f.write(f"{annotation['start']}\\t{annotation['end']}\\t{annotation['score']}\\t")
        f.write(f"{annotation['strand']}\\t{annotation['phase']}\\t{annotation['attributes']}\\n")

# Generate a summary report
profile_counts = defaultdict(int)
class_counts = defaultdict(int)

for key, polyprotein in polyproteins.items():
    profile_counts[polyprotein["profile"]] += 1
    class_counts[polyprotein["completeness_class"]] += 1

# Create DataFrames for the report
profile_df = pd.DataFrame({
    "Profile": list(profile_counts.keys()),
    "Count": list(profile_counts.values())
}).sort_values("Count", ascending=False)

class_df = pd.DataFrame({
    "Completeness Class": list(class_counts.keys()),
    "Count": list(class_counts.values())
}).sort_values("Count", ascending=False)

# Generate HTML report
html = f"""
<html>
<head>
    <title>Polyprotein Annotation Report - {sample_id}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin-bottom: 20px;
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
    <h1>Polyprotein Annotation Report</h1>
    <div class="summary">
        <h2>Summary</h2>
        <p>Sample ID: {sample_id}</p>
        <p>Total polyproteins: {len(polyproteins)}</p>
    </div>
    
    <h2>Polyproteins by Profile</h2>
    {profile_df.to_html(index=False) if not profile_df.empty else "<p>No profiles found.</p>"}
    
    <h2>Polyproteins by Completeness Class</h2>
    {class_df.to_html(index=False) if not class_df.empty else "<p>No completeness classes found.</p>"}
</body>
</html>
"""

with open("${sample_id}_annotation_report.html", "w") as f:
    f.write(html)
EOF
    """
}