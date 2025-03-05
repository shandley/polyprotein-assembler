# CLAUDE.md - Agent Memory Document

## Project: Polyprotein Assembler
This project is a Nextflow workflow for assembling viral polyproteins from metagenomic data.

## Commands

- **Run Pipeline**:
  ```bash
  nextflow run main.nf --reads 'path/to/reads/*_{1,2}.fastq.gz'
  ```

- **Run with Docker**:
  ```bash
  nextflow run main.nf -profile docker --reads 'path/to/reads/*_{1,2}.fastq.gz'
  ```

- **Run with Singularity**:
  ```bash
  nextflow run main.nf -profile singularity --reads 'path/to/reads/*_{1,2}.fastq.gz'
  ```

- **Run with Conda**:
  ```bash
  nextflow run main.nf -profile conda --reads 'path/to/reads/*_{1,2}.fastq.gz'
  ```

- **Build HMM Profiles**:
  ```bash
  hmmbuild profile.hmm aligned_sequences.afa
  cat profile1.hmm profile2.hmm > db/polyprotein_profiles.hmm
  hmmpress db/polyprotein_profiles.hmm
  ```

## Code Style Guidelines

- **Nextflow**:
  - Process names should be ALL_CAPS
  - Parameters should use camelCase
  - Use DSL2 syntax

- **Python**:
  - Follow PEP 8 conventions
  - Use snake_case for functions and variables
  - Use PascalCase for classes
  - Use UPPERCASE for constants
  - Use docstrings for all functions and classes
  - Use specific exceptions in try/except blocks
  - Group imports: standard library, third-party, local
  - Use type hints (Python 3.6+)

- **Shell Scripts**:
  - Use lowercase for variable names
  - Comment complex commands
  - Include error handling with set -e

- **Documentation**:
  - Markdown for all documentation files
  - Include usage examples
  - Document input and output formats