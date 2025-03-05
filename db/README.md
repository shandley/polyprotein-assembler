# HMM Profile Database

This directory should contain the HMM profiles for viral polyproteins.

## Creating HMM Profiles

To create custom HMM profiles for viral polyproteins:

1. Collect polyprotein sequences from reference databases (e.g., RefSeq, GenBank)
2. Filter for viruses of interest (e.g., mammalian viruses, picornaviruses)
3. Cluster sequences to reduce redundancy:
   ```bash
   cd-hit -i all_polyproteins.faa -o clustered_polyproteins.faa -c 0.9
   ```
4. Create multiple sequence alignments:
   ```bash
   mafft --auto clustered_polyproteins.faa > aligned_polyproteins.afa
   ```
5. Build HMM profiles:
   ```bash
   hmmbuild polyprotein_profile.hmm aligned_polyproteins.afa
   ```
6. Combine multiple profiles and prepare for searching:
   ```bash
   cat profile1.hmm profile2.hmm > polyprotein_profiles.hmm
   hmmpress polyprotein_profiles.hmm
   ```

## Example Profiles for Picornavirus Polyproteins

For picornaviruses and related mammalian viruses, consider creating profiles for:

1. Complete polyproteins
2. Conserved domains:
   - 3D RNA-dependent RNA polymerase
   - 2C helicase
   - 3C protease
   - Capsid proteins (VP1-4)
   
## Notes

- Place the final `polyprotein_profiles.hmm` file in this directory
- Ensure the file is indexed using `hmmpress` before running the pipeline
- The pipeline automatically uses this file as the default HMM database