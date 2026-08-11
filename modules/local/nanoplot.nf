#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process nanoplot {
	label "medium"
	publishDir "${params.out_dir}/nanoplot",mode:"copy"
	input:
	tuple val(SampleName),path(fastq)
	output:
	path("${SampleName}_NanoStats_unfilt.txt"),emit:stats_ufilt
	path("${SampleName}_NanoPlot-report_unfilt.html")
	script:
	"""
	NanoPlot --fastq ${fastq} --percentqual -o ${SampleName}_nanoplot_report
	mv ${SampleName}_nanoplot_report/NanoStats.txt ${SampleName}_NanoStats_unfilt.txt
	mv ${SampleName}_nanoplot_report/NanoPlot-report.html ${SampleName}_NanoPlot-report_unfilt.html
	"""
}
