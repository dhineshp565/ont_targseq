#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process minimap2_dehost {
	publishDir "${params.out_dir}/dehost/"
	label "high"
	input:
	tuple val(SampleName),path(fastq)
	path (host_reference)
	output:
	tuple val(SampleName),path ("${SampleName}_host.sam")
	script:
	"""
	minimap2 -ax map-ont -k 15 -w 10 ${host_reference} ${fastq} > ${SampleName}_host.sam
	"""
}

process samtools_dehost {

	publishDir "${params.out_dir}/dehost/"
	label "high"

	input:
	tuple val(SampleName),path (samfile)

	output:
	tuple val(SampleName), path ("${SampleName}_dehosted.fastq.gz"), emit: dehosted_reads
	path("${SampleName}_dehost_stats.txt"), emit: dehost_stats

	script:
	"""
	
	samtools view -b -h ${samfile}|samtools sort > ${SampleName}_hosted.bam
	# Generate dehosting statistics
	samtools stats ${SampleName}_hosted.bam > ${SampleName}_dehost_stats.txt
	# dehost with f 4 flag
	samtools view -b -h -f 4 ${SampleName}_hosted.bam|samtools sort > ${SampleName}_dehosted.bam
	#coverst dehosted bam to fastq
	samtools fastq ${SampleName}_dehosted.bam | gzip > ${SampleName}_dehosted.fastq.gz
	
	
	"""

}

process kraken2_dehost{

	publishDir "${params.out_dir}/dehost/"
	label "high"
	input:
	tuple val(SampleName),path(fastq)
	path (host_db)
	output:
	tuple val(SampleName),path ("${SampleName}_dehosted.fastq.gz"), emit:dehosted_reads
	path ("${SampleName}_kraken_dehost_report.csv"),emit:dehost_stats
	
	script:
	"""
	kraken2 --db ${host_db} --output ${SampleName}_kraken_dehost.csv --report ${SampleName}_kraken_dehost_report.csv --unclassified-out ${SampleName}_dehosted.fastq --classified-out ${SampleName}_bostaurus.fastq.gz ${fastq}
	gzip  ${SampleName}_dehosted.fastq
	"""
}

workflow DEHOST {
	
	take:
	reads           // Channel: tuple val(SampleName), path(fastq)
	host_db  // Path: host genome reference
	
	main:
	// minimap2_dehost(reads, host_reference)
	// samtools_dehost(minimap2_dehost.out)
	kraken2_dehost(reads,host_db)
	
	emit:
	dehosted_reads = kraken2_dehost.out.dehosted_reads
	stats = kraken2_dehost.out.dehost_stats
}