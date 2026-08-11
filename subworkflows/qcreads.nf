#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include local modules
include { make_csv } from '../modules/local/make_csv.nf'
include { merge_filter_fastq } from '../modules/local/merge_filter_fastq.nf'
include { nanoplot } from '../modules/local/nanoplot.nf'
include { porechop } from '../modules/local/porechop.nf'


workflow QCREADS {
    take:
        fastq_path // Channel: input data paths
        qscore // Minimum qscore for filtering
        trim_barcodes // Boolean to trim barcodes
    main:
        data = channel.fromPath(fastq_path)
        
        make_csv(data)
        
        ch_samples = make_csv.out
            .splitCsv(header:true)
            .map { row -> tuple(row.SampleName, row.SamplePath) }

        merge_filter_fastq(ch_samples, qscore)

        //read statistics
        nanoplot(merge_filter_fastq.out.reads)    

        //trim barcodes and adapter sequences
        if (trim_barcodes){
            porechop(merge_filter_fastq.out.reads)
            reads = porechop.out       
        } else {
            reads = merge_filter_fastq.out.reads
        }
    
    emit:
        reads
        read_stats = nanoplot.out.stats_ufilt
        csv = make_csv.out
}