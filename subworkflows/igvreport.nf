#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mapped_ref_bed {
    label "low"
    publishDir "${params.out_dir}/mapped_ref_bed"	
    
    input:
    path (reference)
    
    output:
    path("*.fai")
    path ("*.bed"),emit:bed
    
    script:
    """
    cp ${reference} reference.fasta
    samtools faidx reference.fasta

    awk 'BEGIN {FS=OFS="\t"}; {print \$1,0,\$2}' reference.fasta.fai > reference.bed
    """
}

process bedtools {
    label "high"
    publishDir "${params.out_dir}/bedtools/"
    
    input:
    tuple val(SampleName), path (bam)
    
    output:
    path ("${SampleName}*.bedgraph")
    
    script:
    """
    bedtools genomecov -ibam ${bam} -bga > ${SampleName}.bedgraph
    """
}


process igvreports {
    label "high"
    publishDir "${params.out_dir}/igvreports/", mode: "copy"
    
    input:
    path(reference)
    path(bed)
    path(bedgraph)
    
    output:
    path("*.html")
    
    script:
    """
    create_report ${bed} ${reference} --tracks *.bedgraph --output igv_coverage.html
    """
}

workflow IGVREPORT {
    
    take:
    reference
    bam

    main:
    mapped_ref_bed(reference)
    
    bedtools(bam)
    
    igvreports(reference,mapped_ref_bed.out.bed,bedtools.out.collect())

    emit:
    igv_report = igvreports.out
}