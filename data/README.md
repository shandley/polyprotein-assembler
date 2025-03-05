# Data Directory

Place paired-end metagenomic sequencing data in this directory.

## Expected File Format

The pipeline expects paired-end FASTQ files with the following naming convention:
- `sample_1.fastq.gz` and `sample_2.fastq.gz`
- or `sample_R1.fastq.gz` and `sample_R2.fastq.gz`

## Example Data

For test purposes, you can download example metagenomic datasets from:
- NCBI SRA (https://www.ncbi.nlm.nih.gov/sra)
- MGnify (https://www.ebi.ac.uk/metagenomics/)
- MG-RAST (https://www.mg-rast.org/)

Example command to download test data:
```bash
# Download SRA data using SRA Toolkit
prefetch SRR12345678
fastq-dump --split-files --gzip SRR12345678
mv SRR12345678_1.fastq.gz SRR12345678_2.fastq.gz ./data/
```

## Note

The pipeline will automatically process all paired-end FASTQ files found in this directory that match the specified pattern (default: `*_{1,2}.fastq.gz`).