process POLYPROTEIN_COMPLETION {
    tag "${sample_id}"
    label 'process_low'
    
    publishDir "${params.outdir}/polyprotein_completion", mode: 'copy'
    
    input:
    tuple val(sample_id), path(polyproteins_json)
    
    output:
    tuple val(sample_id), path("${sample_id}_scored_polyproteins.json"), emit: scored_polyproteins
    path "${sample_id}_polyprotein_completeness.png"
    
    script:
    """
    # Score polyproteins based on completeness and other metrics
    python3 << EOF
import json
import os
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

# Load the polyproteins
with open("${polyproteins_json}") as f:
    polyproteins = json.load(f)

# Score each polyprotein
scored_polyproteins = {}
completeness_values = []
profiles = defaultdict(list)

for key, polyprotein in polyproteins.items():
    # Base score is completeness
    base_score = polyprotein["completeness"]
    
    # Add bonus for high fragment count, indicating good support
    fragment_bonus = min(0.1, polyprotein["fragment_count"] * 0.02)
    
    # Penalty for low total score
    score_penalty = 0
    if polyprotein["total_score"] < 100:
        score_penalty = 0.1
    
    # Calculate final score
    final_score = base_score + fragment_bonus - score_penalty
    
    # Classify completeness
    if base_score >= 0.9:
        completeness_class = "Complete"
    elif base_score >= 0.7:
        completeness_class = "Near Complete"
    elif base_score >= 0.5:
        completeness_class = "Partial"
    else:
        completeness_class = "Fragment"
    
    # Store the scored polyprotein
    scored_polyproteins[key] = polyprotein.copy()
    scored_polyproteins[key]["final_score"] = final_score
    scored_polyproteins[key]["completeness_class"] = completeness_class
    
    completeness_values.append(base_score)
    profiles[polyprotein["profile"]].append(base_score)

# Generate a completeness distribution plot
plt.figure(figsize=(10, 6))
plt.hist(completeness_values, bins=20, alpha=0.7, color='blue', edgecolor='black')
plt.xlabel('Completeness')
plt.ylabel('Count')
plt.title('Polyprotein Completeness Distribution')
plt.grid(True, alpha=0.3)
plt.savefig("${sample_id}_polyprotein_completeness.png", dpi=300, bbox_inches='tight')

# Save the scored polyproteins
with open("${sample_id}_scored_polyproteins.json", 'w') as f:
    json.dump(scored_polyproteins, f, indent=2)
EOF
    """
}