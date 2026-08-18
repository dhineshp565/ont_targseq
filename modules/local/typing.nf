#!/usr/bin/env nextflow

process typing {
    publishDir "${params.out_dir}/typing/",mode:"copy"
    label "low"
    
    input:
    tuple val(SampleName),path(consensus)
    path(dbdir)
    
    output:
    path("${SampleName}_abricate.tsv"),emit:abricate
    path("${SampleName}_withseq.tsv"),emit:withseq
    path("${SampleName}_filtered.fasta"),emit:filtered
    
    script:
    """
    typing.py ${SampleName} ${consensus} ${dbdir}
    """
}
