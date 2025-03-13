# Polyprotein Assembler

<div align="center">
  <b>🚧 This project is currently under construction 🚧</b>
</div>
------------------------------------------------------------

A Nextflow workflow to assemble and reconstruct viral polyproteins from metagenomic data, with a focus on mammalian viruses like picornaviruses.

## Overview

This workflow is designed to identify, assemble, and reconstruct viral polyproteins from metagenomic sequencing data. The pipeline employs a custom HMM profile-based approach for sensitive detection of viral polyprotein fragments, followed by targeted reassembly and fragment stitching to produce more complete polyproteins.

## Key Features

- Quality control and pre-processing of metagenomic reads
- Metagenomic assembly with flexible assembly tool selection
- Sensitive detection of viral polyproteins using custom HMM profiles
- Targeted reassembly of contigs containing polyprotein fragments
- Fragment stitching to reconstruct more complete polyproteins
- Scoring and classification of polyproteins by completeness
- Functional annotation and visualization

## Workflow Steps

1. **Read QC**: Quality control and adapter trimming using FastQC and Trimmomatic
2. **Assembly**: Metagenomic assembly using MEGAHIT or SPAdes
3. **Translation**: 6-frame translation of assembled contigs
4. **HMM Search**: Identification of polyprotein regions using custom HMM profiles
5. **Targeted Reassembly**: Improved assembly around polyprotein-containing regions
6. **Fragment Stitching**: Reconstruction of polyproteins by connecting fragments
7. **Polyprotein Completion**: Scoring and classification of polyproteins by completeness
8. **Annotation**: Functional annotation of reconstructed polyproteins

## Installation

```bash
# Clone the repository
git clone https://github.com/scotthandley/polyprotein-assembler.git
cd polyprotein-assembler

# Create required directories
mkdir -p db assets data
```

## Usage

```bash
# Basic usage
nextflow run main.nf --reads 'path/to/reads/*_{1,2}.fastq.gz'

# Specify output directory
nextflow run main.nf --reads 'path/to/reads/*_{1,2}.fastq.gz' --outdir my_results

# Use SPAdes for assembly
nextflow run main.nf --reads 'path/to/reads/*_{1,2}.fastq.gz' --assembly_tool spades

# Run with Docker
nextflow run main.nf -profile docker --reads 'path/to/reads/*_{1,2}.fastq.gz'

# Run with Conda/Mamba (recommended for all systems)
nextflow run main.nf -profile conda --reads 'path/to/reads/*_{1,2}.fastq.gz'

# For Apple Silicon (M1/M2/M3) optimized performance
nextflow run main.nf -profile apple_silicon --reads 'path/to/reads/*_{1,2}.fastq.gz'
```

### Apple Silicon Optimization

The pipeline includes specific optimizations for Apple Silicon (M1/M2/M3) Macs:

1. **Installation on Apple Silicon**:
   ```bash
   # Create conda environment optimized for Apple Silicon
   CONDA_SUBDIR=osx-arm64 mamba env create -f conda/conda-apple-silicon.yml
   ```

2. **Performance Optimizations**:
   - Use of `pyhmmer` for native ARM64 HMMER searches (~2-4x faster than emulated HMMER)
   - Optimized memory usage for memory-intensive steps (assembly, HMM search)
   - Native ARM64 binaries for all tools when using the `apple_silicon` profile

3. **Running with Apple Silicon Optimizations**:
   ```bash
   # Use the apple_silicon profile for best performance
   nextflow run main.nf -profile apple_silicon --reads 'path/to/reads/*_{1,2}.fastq.gz'
   ```

## Custom HMM Profiles

To use the pipeline with custom HMM profiles:

1. Create HMM profiles from aligned polyprotein sequences using HMMER:
   ```bash
   hmmbuild my_profile.hmm aligned_polyproteins.afa
   ```

2. Combine multiple profiles into a single database:
   ```bash
   cat profile1.hmm profile2.hmm > db/polyprotein_profiles.hmm
   hmmpress db/polyprotein_profiles.hmm
   ```

3. Run the pipeline with your custom profiles:
   ```bash
   nextflow run main.nf --reads 'path/to/reads/*_{1,2}.fastq.gz' --hmm_db /path/to/custom_profiles.hmm
   ```

## Output

The pipeline generates the following outputs in the results directory:

- `fastqc/`: Read quality reports
- `trimmed/`: Trimmed reads
- `assembly/`: Assembled contigs
- `hmm_search/`: HMM search results
- `targeted_reassembly/`: Reassembled contigs around polyprotein regions
- `fragment_stitching/`: Fragment stitching results and reports
- `polyprotein_completion/`: Polyprotein completeness scores and reports
- `annotation/`: Annotated polyproteins and functional reports

## Requirements

- Nextflow (>=21.04.0)
- Java (>=11)
- Docker, Conda/Mamba, or Singularity (optional, for containerized execution)

The following tools are required (automatically installed with Docker/Conda profiles):
- FastQC
- Trimmomatic
- MEGAHIT/SPAdes
- EMBOSS (transeq)
- HMMER (or pyhmmer on Apple Silicon)
- Python (>=3.9) with Biopython, pandas, matplotlib, and networkx
- BWA
- Samtools
- Seqtk
- CD-HIT (for clustering sequences when building HMM profiles)
- MAFFT (for creating alignments when building HMM profiles)

### System Requirements

- **Recommended Memory**: 16GB+ RAM (32GB+ for large metagenomes)
- **Disk Space**: 100GB+ free space for temporary and output files
- **CPU**: 8+ cores recommended
- **Supported Platforms**:
  - Linux (x86_64)
  - macOS (Intel and Apple Silicon)
  - Windows (via WSL2 or Docker)
  
For Apple Silicon (M1/M2/M3) Macs, we recommend using the `apple_silicon` profile for optimal performance.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
