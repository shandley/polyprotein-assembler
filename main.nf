#!/usr/bin/env nextflow

/*
========================================================================================
    Polyprotein Assembler Pipeline
========================================================================================
    A Nextflow workflow for assembling viral polyproteins from metagenomic data
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Print pipeline header
log.info """
===========================================
 Polyprotein Assembler v1.0
===========================================
"""

// Define parameters
params.reads = "$projectDir/data/*_{1,2}.fastq.gz"
params.outdir = "results"
params.threads = 8
params.hmm_db = "$projectDir/db/polyprotein_profiles.hmm"
params.assembly_tool = "megahit"
params.min_contig_length = 500
params.help = false

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Include processes
include { FASTQC } from './modules/fastqc'
include { TRIMMING } from './modules/trimming'
include { ASSEMBLY } from './modules/assembly'
include { TRANSLATE_CONTIGS } from './modules/translate'
include { HMM_SEARCH } from './modules/hmm_search'
include { POLYPROTEIN_COMPLETION } from './modules/polyprotein_completion'
include { TARGETED_REASSEMBLY } from './modules/targeted_reassembly'
include { FRAGMENT_STITCHING } from './modules/fragment_stitching'
include { ANNOTATION } from './modules/annotation'

// Define input channels
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

// Main workflow
workflow {
    FASTQC(read_pairs_ch)
    TRIMMING(read_pairs_ch)
    ASSEMBLY(TRIMMING.out.trimmed_reads)
    TRANSLATE_CONTIGS(ASSEMBLY.out.contigs)
    HMM_SEARCH(TRANSLATE_CONTIGS.out.translated_contigs, params.hmm_db)
    
    // Separate full and partial polyprotein hits
    HMM_SEARCH.out.hmm_hits
        .map { it -> [it[0], it[1].findAll { hit -> hit.score >= params.full_polyprotein_threshold }] }
        .set { full_polyproteins_ch }
    
    HMM_SEARCH.out.hmm_hits
        .map { it -> [it[0], it[1].findAll { hit -> hit.score < params.full_polyprotein_threshold && hit.score >= params.partial_polyprotein_threshold }] }
        .set { partial_polyproteins_ch }
    
    // Process partial polyproteins
    TARGETED_REASSEMBLY(partial_polyproteins_ch, TRIMMING.out.trimmed_reads)
    FRAGMENT_STITCHING(
        TARGETED_REASSEMBLY.out.reassembled_contigs,
        full_polyproteins_ch
    )
    
    // Combine full and improved partial polyproteins
    full_polyproteins_ch
        .mix(FRAGMENT_STITCHING.out.improved_polyproteins)
        .set { combined_polyproteins_ch }
    
    // Score and annotate polyproteins
    POLYPROTEIN_COMPLETION(combined_polyproteins_ch)
    ANNOTATION(POLYPROTEIN_COMPLETION.out.scored_polyproteins)
}

def helpMessage() {
    log.info"""
    =========================================
    Polyprotein Assembler Pipeline
    =========================================
    Usage:
    nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --outdir results
    
    Required arguments:
      --reads                Path to input data (must be quoted)
    
    Optional arguments:
      --outdir               Output directory (default: 'results')
      --threads              Number of CPUs to use (default: 8)
      --hmm_db               Path to HMM database (default: 'db/polyprotein_profiles.hmm')
      --assembly_tool        Assembly tool to use: 'megahit' or 'spades' (default: 'megahit')
      --min_contig_length    Minimum contig length to retain (default: 500)
    """
}