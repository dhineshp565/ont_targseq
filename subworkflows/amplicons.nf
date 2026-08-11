#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { minimap2 } from '../modules/local/minimap2.nf'
include {splitbam } from '../modules/local/splitbam.nf'
include { medaka } from '../modules/local/medaka.nf'

workflow AMPLICONS {

    take:
    reads
    reference
    primerbed
    readcount
    consensus_mode
    qscore

    main:
        minimap2(reference, reads)
        splitbam(minimap2.out, primerbed, readcount, consensus_mode, qscore)
         // Join reads and consensus for medaka
        medaka_input=   splitbam.out.consensus
                        .join(reads)
                        .map { sample, consensus, fastq -> tuple(sample, fastq, consensus) }
        
        medaka(medaka_input)

    emit:
    consensus   = medaka.out.consensus
    cons_only   = splitbam.out.cons_only
    stats       = splitbam.out.unfilt_stats
    idxstats    = splitbam.out.idxstats
    mapped      = splitbam.out.mapped
    bam         = splitbam.out.target_bam
   



}